%% ========================================================================
%  radar_beam_sim.m
%  송신빔 1개 + 수신빔 9개(중앙 1 + 외곽 8, 3x3 사각배치) 레이더 빔 시뮬레이션
% -------------------------------------------------------------------------
%  기능
%   - 각도공간(az/el, deg)에서 안테나 크기/파장 기반 sinc 빔패턴 모의
%   - 표적(위치 az/el, RCS) 모의
%   - 송신이득 x 각 수신빔 이득(2-way) x RCS + 잡음 으로 빔별 신호레벨 생성
%   - 추정 알고리즘 3종 : max-amp / centroiding / 가로,세로 parabolic peak
%   - 단일표적 몬테카를로 추정오차 분석 (RMSE, 산점도, 히스토그램,
%                                       SNR 대비/오프셋 대비 RMSE)
%   - 다중표적 운용 시 추정오차 분석 (표적개수, 표적간격 대비 RMSE)
%
%  실행 : >> radar_beam_sim
%  의존성 : 없음 (sinc 는 로컬함수로 직접 구현, 툴박스 불필요)
%
%  ※ 빔 위치/간격은 아래 cfg 와 make_rx_beams() 에서 자유롭게 수정 가능.
% =========================================================================
function radar_beam_sim
clc; close all; rng(2026);   % 재현성

%% ------------------------- 1) 설정(Config) ------------------------------
cfg.c       = 3e8;            % 광속 [m/s]
cfg.f       = 10e9;          % 주파수 [Hz] (X-band 예시)
cfg.lambda  = cfg.c/cfg.f;   % 파장 [m]

% 안테나 개구면 크기 [m] (가로=az, 세로=el). sinc 패턴의 빔폭을 결정.
cfg.Dtx_az  = 0.30;          % 송신 안테나 가로 크기
cfg.Dtx_el  = 0.30;          % 송신 안테나 세로 크기
cfg.Drx_az  = 0.30;          % 수신 안테나 가로 크기 (각 수신빔 동일 가정)
cfg.Drx_el  = 0.30;          % 수신 안테나 세로 크기

% 수신빔 배치 : 중앙 1 + 외곽 8 (3x3). 간격 조정 가능.
cfg.beam_spacing = 1.0;      % 빔 간격 [deg]  (외곽빔 오프셋)
cfg.tx_point     = [0 0];    % 송신빔 지향각 [az el] (deg)

% 기준 SNR : 중앙빔 보어사이트에 RCS=RCS_ref 표적이 있을 때의 SNR [dB]
cfg.SNR0_dB   = 30;
cfg.RCS_ref   = 1.0;         % 기준 RCS [m^2]
cfg.RCS_def   = 1.0;         % 단일 시나리오 기본 표적 RCS [m^2]

% 몬테카를로 설정
cfg.MC_N      = 3000;        % 시행 횟수
cfg.fov       = cfg.beam_spacing;   % 표적 랜덤 발생 범위 ±fov [deg] (모호성 영역내)

% 수신빔 생성 (위치/인덱스)
[beam, idx] = make_rx_beams(cfg.beam_spacing);

% 빔폭(참고) 출력 : sinc 개구면 -3dB 빔폭 ≈ 0.886*lambda/D
bw_tx = 0.886*cfg.lambda/cfg.Dtx_az*180/pi;
bw_rx = 0.886*cfg.lambda/cfg.Drx_az*180/pi;
fprintf('파장=%.4f m | 송신 -3dB 빔폭 ~%.3f deg | 수신 -3dB 빔폭 ~%.3f deg | 빔간격=%.2f deg\n',...
    cfg.lambda, bw_tx, bw_rx, cfg.beam_spacing);
fprintf('수신빔 %d개 (중앙1 + 외곽8)\n\n', size(beam,1));

%% ----------------- 2) 단일표적 예시 시나리오 + 시각화 -------------------
demo_single_scenario(cfg, beam, idx);

%% ----------------- 3) 단일표적 몬테카를로 오차분석 ----------------------
mc_single_target(cfg, beam, idx);

%% ----------------- 4) 다중표적 운용 오차분석 ---------------------------
mc_multi_target(cfg, beam, idx);

fprintf('\n완료. 생성된 Figure 들을 확인하세요.\n');
end


%% ========================================================================
%  수신빔 배치 생성 : 중앙 1 + 외곽 8 (3x3 사각배치)
%  반환 : beam(Nx2) = [az el] (deg),  idx 구조체(center/E/W/N/S)
%  ── 빔 위치를 임의로 바꾸고 싶으면 아래 beam 행렬을 직접 수정.
% ========================================================================
function [beam, idx] = make_rx_beams(s)
% 3x3 격자 : 위(+el)부터 아래, 좌(-az)에서 우 순서
az = [-s 0 s];
el = [ s 0 -s];
beam = zeros(9,2); k = 0;
for ie = 1:3
    for ia = 1:3
        k = k+1;
        beam(k,:) = [az(ia) el(ie)];
    end
end
% 추정에 쓰는 십자(cross) 인덱스 찾기 (parabolic peak 용)
idx.C = find_beam(beam, [ 0  0]);   % 중앙
idx.E = find_beam(beam, [ s  0]);   % 동(+az)
idx.W = find_beam(beam, [-s  0]);   % 서(-az)
idx.N = find_beam(beam, [ 0  s]);   % 북(+el)
idx.S = find_beam(beam, [ 0 -s]);   % 남(-el)
idx.spacing = s;
end

function k = find_beam(beam, p)
[~,k] = min(sum((beam - p).^2,2));
end


%% ========================================================================
%  sinc 기반 직사각 개구면 안테나 전력패턴 (one-way, 선형, 피크=1)
%  daz,del : 빔 보어사이트 기준 각도오차 [deg]
% ========================================================================
function G = antenna_pattern(daz, del, Daz, Del, lambda)
ua = (Daz/lambda).*sind(daz);
ue = (Del/lambda).*sind(del);
G  = (mysinc(ua).*mysinc(ue)).^2;   % 전력 = (전압패턴)^2
end

function y = mysinc(x)
% 정규화 sinc = sin(pi x)/(pi x), x=0 -> 1
y = ones(size(x));
n = (x~=0);
y(n) = sin(pi*x(n))./(pi*x(n));
end


%% ========================================================================
%  빔별 측정 신호 생성
%  targets : K x 3 = [az el RCS]
%  반환 amp : N x 1 (각 수신빔의 측정 진폭, 잡음포함)
% ========================================================================
function amp = gen_beam_signals(targets, cfg, beam)
N = size(beam,1);
Vsig = zeros(N,1);

% 기준전력 Pref=1 기준 잡음전력
Pref       = 1.0;
noisePow   = Pref / 10^(cfg.SNR0_dB/10);
noiseSigma = sqrt(noisePow/2);   % I,Q 각 채널 표준편차

for t = 1:size(targets,1)
    tAz = targets(t,1); tEl = targets(t,2); rcs = targets(t,3);
    phi = 2*pi*rand;             % 표적별 임의 위상 (독립 산란체)

    % 송신 이득 (송신빔 지향각 기준)
    Gtx = antenna_pattern(tAz-cfg.tx_point(1), tEl-cfg.tx_point(2), ...
                          cfg.Dtx_az, cfg.Dtx_el, cfg.lambda);

    for b = 1:N
        Grx = antenna_pattern(tAz-beam(b,1), tEl-beam(b,2), ...
                              cfg.Drx_az, cfg.Drx_el, cfg.lambda);
        Prx = Pref * Gtx * Grx * (rcs/cfg.RCS_ref);   % 2-way 상대 수신전력
        Vsig(b) = Vsig(b) + sqrt(Prx)*exp(1i*phi);    % 전압 합산
    end
end

% 열잡음 추가 후 진폭(포락선) 측정
noise = noiseSigma*(randn(N,1)+1i*randn(N,1));
amp   = abs(Vsig + noise);
end


%% ========================================================================
%  추정 알고리즘 3종.  est = [az el] (deg)
% ========================================================================
% (1) Max amplitude : 최대 진폭 빔의 위치
function est = estimate_maxamp(amp, beam)
[~,k] = max(amp);
est = beam(k,:);
end

% (2) Centroiding : 전력가중 무게중심 (전 빔)
function est = estimate_centroid(amp, beam)
p = amp.^2;                     % 전력 가중
% 잡음바닥 영향 완화를 위해 최대전력 대비 임계치 이하 제거
p = max(p - 0.05*max(p), 0);
if sum(p)==0, est = [0 0]; return; end
est = [sum(p.*beam(:,1)) sum(p.*beam(:,2))] / sum(p);
end

% (3) Peak estimation : 가로(W,C,E)/세로(S,C,N) 3점 포물선 보간
function est = estimate_peak(amp, beam, idx)
p  = amp.^2;
s  = idx.spacing;
% 가로(az)
daz = parabolic_offset(p(idx.W), p(idx.C), p(idx.E));
% 세로(el)
del = parabolic_offset(p(idx.S), p(idx.C), p(idx.N));
estAz = beam(idx.C,1) + daz*s;   % S->C->N 순(=-el->+el)이므로 동일부호
estEl = beam(idx.C,2) + del*s;
est = [estAz estEl];
end

function d = parabolic_offset(yL, yC, yR)
% 3점 포물선 정점 오프셋 (-1..1), +방향(R)쪽이 양수
den = (yL - 2*yC + yR);
if abs(den) < eps
    d = 0;
else
    d = 0.5*(yL - yR)/den;
end
d = max(min(d,1),-1);            % 격자 밖으로 발산 방지
end


%% ========================================================================
%  단일표적 예시 시나리오 + 빔패턴/진폭/추정결과 시각화
% ========================================================================
function demo_single_scenario(cfg, beam, idx)
tgt = [0.35, -0.25, cfg.RCS_def];     % 예시 표적 [az el RCS]
amp = gen_beam_signals(tgt, cfg, beam);

eM = estimate_maxamp(amp, beam);
eC = estimate_centroid(amp, beam);
eP = estimate_peak(amp, beam, idx);

fprintf('[예시 시나리오] 참 표적=(%.3f, %.3f) deg\n', tgt(1), tgt(2));
fprintf('   max-amp  : (%.3f, %.3f)  err=%.3f deg\n', eM, norm(eM-tgt(1:2)));
fprintf('   centroid : (%.3f, %.3f)  err=%.3f deg\n', eC, norm(eC-tgt(1:2)));
fprintf('   peak     : (%.3f, %.3f)  err=%.3f deg\n\n', eP, norm(eP-tgt(1:2)));

% --- Fig1 : 빔패턴(-3dB 윤곽) + 빔중심 + 표적 ---
figure('Name','Fig1 빔배치/빔패턴','Color','w');
g = linspace(-2*cfg.beam_spacing, 2*cfg.beam_spacing, 241);
[AZ,EL] = meshgrid(g,g);
hold on;
% 송신빔 -3dB
Gt = antenna_pattern(AZ-cfg.tx_point(1), EL-cfg.tx_point(2), cfg.Dtx_az,cfg.Dtx_el,cfg.lambda);
contour(AZ,EL,10*log10(Gt+eps),[-3 -3],'k--','LineWidth',1.5);
% 각 수신빔 -3dB
for b = 1:size(beam,1)
    Gr = antenna_pattern(AZ-beam(b,1), EL-beam(b,2), cfg.Drx_az,cfg.Drx_el,cfg.lambda);
    contour(AZ,EL,10*log10(Gr+eps),[-3 -3],'b','LineWidth',1.0);
end
plot(beam(:,1),beam(:,2),'bo','MarkerFaceColor','b','MarkerSize',5);
plot(tgt(1),tgt(2),'rp','MarkerFaceColor','r','MarkerSize',16);
text(beam(idx.C,1),beam(idx.C,2),'  C','Color','b');
axis equal; grid on; xlabel('Azimuth [deg]'); ylabel('Elevation [deg]');
title('수신빔 배치(파랑 -3dB) / 송신빔(검정점선) / 표적(별)');
legend('Tx -3dB','Rx -3dB','location','northeastoutside');

% --- Fig2 : 빔별 측정진폭(빔위치에 색/크기 표시) + 추정결과 ---
figure('Name','Fig2 빔진폭/추정','Color','w');
sc = 60 + 600*(amp/max(amp));        % 진폭 비례 마커크기
scatter(beam(:,1),beam(:,2),sc,amp,'filled'); hold on;
colorbar; axis equal;
xlim([-1.6 1.6]*cfg.beam_spacing); ylim([-1.6 1.6]*cfg.beam_spacing);
plot(tgt(1),tgt(2),'rp','MarkerFaceColor','r','MarkerSize',16);
plot(eM(1),eM(2),'ws','MarkerSize',12,'LineWidth',2);
plot(eC(1),eC(2),'g^','MarkerSize',11,'LineWidth',2);
plot(eP(1),eP(2),'mo','MarkerSize',11,'LineWidth',2);
xlabel('Azimuth [deg]'); ylabel('Elevation [deg]');
title('빔별 측정진폭 + 추정 (별:참, □:max, △:centroid, ○:peak)');
end


%% ========================================================================
%  단일표적 몬테카를로 오차분석
% ========================================================================
function mc_single_target(cfg, beam, idx)
N = cfg.MC_N; f = cfg.fov;
true_pos = (2*rand(N,2)-1)*f;            % ±fov 내 균일분포
err = struct('max',zeros(N,2),'cen',zeros(N,2),'pk',zeros(N,2));

for n = 1:N
    tgt = [true_pos(n,:) cfg.RCS_def];
    amp = gen_beam_signals(tgt, cfg, beam);
    err.max(n,:) = estimate_maxamp(amp,beam)      - true_pos(n,:);
    err.cen(n,:) = estimate_centroid(amp,beam)    - true_pos(n,:);
    err.pk (n,:) = estimate_peak(amp,beam,idx)    - true_pos(n,:);
end

rmse = @(e) sqrt(mean(sum(e.^2,2)));     % 2D 합성 RMSE
fprintf('[단일표적 MC, N=%d, SNR0=%ddB] 합성 RMSE [deg]\n',N,cfg.SNR0_dB);
fprintf('   max-amp =%.4f | centroid=%.4f | peak=%.4f\n\n',...
    rmse(err.max), rmse(err.cen), rmse(err.pk));

% --- Fig3 : 오차 산점도 ---
figure('Name','Fig3 오차 산점도','Color','w');
plot(err.max(:,1),err.max(:,2),'.','Color',[.85 .3 .3]); hold on;
plot(err.cen(:,1),err.cen(:,2),'.','Color',[.2 .6 .2]);
plot(err.pk (:,1),err.pk (:,2),'.','Color',[.2 .2 .8]);
axis equal; grid on; xlabel('Az error [deg]'); ylabel('El error [deg]');
title('단일표적 추정오차 산점도'); legend('max-amp','centroid','peak');

% --- Fig4 : 오차 히스토그램(합성거리) ---
figure('Name','Fig4 오차 히스토그램','Color','w');
em = sqrt(sum(err.max.^2,2)); ec = sqrt(sum(err.cen.^2,2)); ep = sqrt(sum(err.pk.^2,2));
edges = linspace(0, max([em;ec;ep]), 40);
histogram(em,edges,'FaceAlpha',.4); hold on;
histogram(ec,edges,'FaceAlpha',.4);
histogram(ep,edges,'FaceAlpha',.4);
grid on; xlabel('|오차| [deg]'); ylabel('빈도');
title('추정오차 크기 분포'); legend('max-amp','centroid','peak');

% --- Fig5 : SNR 대비 RMSE ---
snr_list = 0:5:40;
R = zeros(numel(snr_list),3);
cfg2 = cfg; cfg2.MC_N = max(800, round(cfg.MC_N/3));
for is = 1:numel(snr_list)
    cfg2.SNR0_dB = snr_list(is);
    R(is,:) = quick_rmse(cfg2, beam, idx);
end
figure('Name','Fig5 SNR 대비 RMSE','Color','w');
plot(snr_list,R(:,1),'-o', snr_list,R(:,2),'-s', snr_list,R(:,3),'-^','LineWidth',1.5);
grid on; xlabel('기준 SNR [dB]'); ylabel('합성 RMSE [deg]');
title('SNR 대비 추정 RMSE'); legend('max-amp','centroid','peak');

% --- Fig6 : 보어사이트 오프셋 대비 RMSE ---
roff = sqrt(sum(true_pos.^2,2));
nb = 8; be = linspace(0, f*sqrt(2), nb+1); bc = (be(1:end-1)+be(2:end))/2;
Rb = nan(nb,3);
for k = 1:nb
    m = roff>=be(k) & roff<be(k+1);
    if nnz(m) > 10
        Rb(k,1)=rmse(err.max(m,:)); Rb(k,2)=rmse(err.cen(m,:)); Rb(k,3)=rmse(err.pk(m,:));
    end
end
figure('Name','Fig6 오프셋 대비 RMSE','Color','w');
plot(bc,Rb(:,1),'-o', bc,Rb(:,2),'-s', bc,Rb(:,3),'-^','LineWidth',1.5);
grid on; xlabel('표적-보어사이트 오프셋 [deg]'); ylabel('합성 RMSE [deg]');
title('표적 위치(오프셋)별 추정 RMSE'); legend('max-amp','centroid','peak');
end

function R = quick_rmse(cfg, beam, idx)
N = cfg.MC_N; f = cfg.fov;
tp = (2*rand(N,2)-1)*f; e = zeros(N,3,2);
for n=1:N
    tgt=[tp(n,:) cfg.RCS_def]; amp=gen_beam_signals(tgt,cfg,beam);
    e(n,1,:)=estimate_maxamp(amp,beam)-tp(n,:);
    e(n,2,:)=estimate_centroid(amp,beam)-tp(n,:);
    e(n,3,:)=estimate_peak(amp,beam,idx)-tp(n,:);
end
R = squeeze(sqrt(mean(sum(e.^2,3),1)))';
end


%% ========================================================================
%  다중표적 운용 오차분석
%   (주의) 단일 빔클러스터(9빔)로는 한 클러스터 내 다중표적의 완전분리가
%   원리적으로 제한됨. 여기서는 (a)표적개수, (b)표적간격 증가에 따른
%   '최강표적' 추정오차 열화를 정량화.
% ========================================================================
function mc_multi_target(cfg, beam, idx)
f = cfg.fov; N = max(800, round(cfg.MC_N/2));

% --- (A) 표적 개수 대비 RMSE (최강표적 기준) ---
Klist = 1:5;
RA = zeros(numel(Klist),3);
for ik = 1:numel(Klist)
    K = Klist(ik); e = zeros(N,3,2);
    for n=1:N
        T = [ (2*rand(K,2)-1)*f , cfg.RCS_def*(0.5+rand(K,1)) ]; % 랜덤 RCS
        amp = gen_beam_signals(T, cfg, beam);
        [~,ks] = max(T(:,3)); ref = T(ks,1:2);  % 최강표적을 평가기준으로
        e(n,1,:)=estimate_maxamp(amp,beam)-ref;
        e(n,2,:)=estimate_centroid(amp,beam)-ref;
        e(n,3,:)=estimate_peak(amp,beam,idx)-ref;
    end
    RA(ik,:)=squeeze(sqrt(mean(sum(e.^2,3),1)))';
end
figure('Name','Fig7 표적개수 대비 RMSE','Color','w');
plot(Klist,RA(:,1),'-o',Klist,RA(:,2),'-s',Klist,RA(:,3),'-^','LineWidth',1.5);
grid on; xlabel('동시 표적 개수'); ylabel('최강표적 합성 RMSE [deg]');
title('다중표적: 표적개수 대비 추정오차'); legend('max-amp','centroid','peak');

% --- (B) 2표적 간격 대비 RMSE (동일 RCS) ---
sep = linspace(0.1, 2.5*cfg.beam_spacing, 12);
RB = zeros(numel(sep),3);
for is = 1:numel(sep)
    d = sep(is); e = zeros(N,3,2);
    for n=1:N
        c0  = (2*rand(1,2)-1)*f*0.4;        % 쌍 중심
        ang = 2*pi*rand;                    % 임의 방향
        off = d/2*[cos(ang) sin(ang)];
        t1 = c0+off; t2 = c0-off;           % 평가기준=t1
        T  = [t1 cfg.RCS_def; t2 cfg.RCS_def];
        amp = gen_beam_signals(T, cfg, beam);
        e(n,1,:)=estimate_maxamp(amp,beam)-t1;
        e(n,2,:)=estimate_centroid(amp,beam)-t1;
        e(n,3,:)=estimate_peak(amp,beam,idx)-t1;
    end
    RB(is,:)=squeeze(sqrt(mean(sum(e.^2,3),1)))';
end
figure('Name','Fig8 2표적 간격 대비 RMSE','Color','w');
plot(sep,RB(:,1),'-o',sep,RB(:,2),'-s',sep,RB(:,3),'-^','LineWidth',1.5);
grid on; xlabel('2표적 각도간격 [deg]'); ylabel('표적1 합성 RMSE [deg]');
title('다중표적: 2표적 간격 대비 추정오차'); legend('max-amp','centroid','peak');

fprintf('[다중표적 MC] 표적개수1->5 RMSE(max/cen/peak):\n');
disp(array2table(RA,'VariableNames',{'maxamp','centroid','peak'},...
     'RowNames',compose('K=%d',Klist)));
end
