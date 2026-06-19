# 빔 기반 표적 위치추정 알고리즘 이론 설명

본 문서는 `radar_beam_sim.m` 에 구현된 3가지 표적 각도위치 추정 알고리즘
(**Max-Amplitude**, **Centroiding**, **Parabolic Peak Estimation**)의
이론적 배경을 정리한다.

---

## 0. 공통 신호 모델

### 0.1 안테나 빔패턴 (sinc 개구면 모델)

크기 `D`(가로 `D_az`, 세로 `D_el`)인 균일 직사각 개구면(uniform rectangular
aperture)의 원거리장(far-field) 전압 패턴은 개구면 전류분포의 푸리에 변환으로
주어지며, 정규화하면 sinc 형태가 된다. **단방향(one-way) 전력 패턴**은

$$
G(\theta_{az},\theta_{el}) =
\Big[\operatorname{sinc}\!\Big(\tfrac{D_{az}}{\lambda}\sin\theta_{az}\Big)\Big]^2
\Big[\operatorname{sinc}\!\Big(\tfrac{D_{el}}{\lambda}\sin\theta_{el}\Big)\Big]^2,
\qquad
\operatorname{sinc}(x)=\frac{\sin(\pi x)}{\pi x}.
$$

- 첫 번째 null: $\dfrac{D}{\lambda}\sin\theta = 1 \Rightarrow
  \theta_{\text{null}} \approx \sin^{-1}(\lambda/D)$
- 반전력 빔폭(HPBW): $\theta_{3\text{dB}} \approx 0.886\,\dfrac{\lambda}{D}$ [rad]

코드에서 `antenna_pattern()` 이 이 식을 구현한다 (전압패턴의 제곱 = 전력패턴).

### 0.2 빔 배치

송신빔 1개(보어사이트 지향)와 수신빔 9개(중앙 1 + 외곽 8, $3\times3$ 격자,
간격 $s$)를 사용한다. 빔 $b$의 지향각을 $(\phi_b,\psi_b)$ 라 한다.

### 0.3 빔별 수신 신호 (2-way)

표적이 $(\theta_{az},\theta_{el})$, RCS $\sigma$ 일 때 빔 $b$의 수신전력은
송신이득 × 수신이득 × RCS 의 2-way 곱으로 모델링한다.

$$
P_b = P_0\, G_{tx}(\theta_{az},\theta_{el})\,
        G_{rx}(\theta_{az}-\phi_b,\ \theta_{el}-\psi_b)\,
        \frac{\sigma}{\sigma_{ref}}
$$

측정 진폭(포락선)은 복소 열잡음 $n_b \sim \mathcal{CN}(0,\sigma_n^2)$ 을 더해

$$
A_b = \big|\sqrt{P_b}\,e^{j\varphi} + n_b\big|,
\qquad
\text{SNR}_0 = \frac{P_0}{\sigma_n^2}\ \ (\text{중앙빔 보어사이트 기준}).
$$

추정기들은 빔별 측정값 $\{A_b\}$ (또는 전력 $A_b^2$)만으로 표적의
$(\theta_{az},\theta_{el})$ 를 추정한다.

---

## 1. Max-Amplitude (최대 진폭)

### 원리
가장 큰 신호를 받은 빔의 지향각을 표적 위치로 추정한다.

$$
\hat b = \arg\max_b A_b,
\qquad
(\hat\theta_{az},\hat\theta_{el}) = (\phi_{\hat b},\psi_{\hat b}).
$$

### 특성
- **장점**: 가장 단순, 계산량 최소, 잡음/간섭에 강건(robust), 모호성 없음.
- **단점**: 추정 해상도가 **빔 격자 간격 $s$ 로 양자화**됨. 추정 오차는
  격자 셀 내 균일분포에 가까워, 1차원 RMSE는 대략

$$
\text{RMSE}_{1D} \approx \frac{s}{\sqrt{12}}
$$

  (셀 폭 $s$ 인 균일분포의 표준편차). 2차원 합성은 $\sqrt2$ 배.
- 빔폭 대비 격자가 촘촘할수록 정확하지만 빔 수가 늘어난다.
- 본질적으로 **편향(bias)** 은 없으나 **분해능(granularity) 한계**가 지배적.

### 용도
초기 탐지/조대(coarse) 추정, 또는 다른 정밀 추정기의 초기 빔 선택.

---

## 2. Centroiding (전력 무게중심)

### 원리
빔 진폭(전력)을 가중치로 한 빔 지향각의 무게중심을 표적 위치로 본다.
여러 빔에 걸친 표적 에너지의 "질량중심"을 구하는 개념이다.

$$
\hat\theta_{az} = \frac{\sum_b w_b\,\phi_b}{\sum_b w_b},
\qquad
\hat\theta_{el} = \frac{\sum_b w_b\,\psi_b}{\sum_b w_b},
\qquad w_b = A_b^2 .
$$

코드에서는 잡음바닥 영향을 줄이기 위해
$w_b = \max(A_b^2 - \kappa\max_b A_b^2,\,0)$ 의 임계(threshold) 처리를
적용한다($\kappa \approx 0.05$). 광학 영상의 무게중심(centroid)·
도심(centroiding) 기법과 동일한 사상이다.

### 특성
- **장점**: 모든 빔 정보를 사용 → 연속적인(서브-격자) 추정값,
  잡음 평균화 효과로 SNR이 높을 때 max-amp보다 정밀.
- **단점 — 계통오차(bias)**: 빔 응답이 좌우 대칭이 아닌 영역(표적이
  격자 가장자리·보어사이트에서 멀어질 때)에서 추정값이 배열 중심 쪽으로
  **끌리는 편향**이 생긴다. 유한 개 빔의 가중합이므로 시야(FOV) 가장자리에서
  특히 두드러진다(코드의 Fig6: 오프셋 대비 RMSE 증가).
- **가중치 선택**: 전력($A^2$) 가중이 일반적이나, 임계치/가중지수에 따라
  편향-분산 절충(bias–variance trade-off)이 달라진다.

### 용도
SNR이 충분하고 표적이 빔 클러스터 중앙부에 있을 때의 정밀 각도추정.
임계처리/배경제거가 정확도에 중요.

---

## 3. Parabolic Peak Estimation (가로/세로 3점 포물선 보간)

### 원리
빔 응답의 주엽(main lobe)은 피크 근방에서 2차(포물선) 함수로 근사된다.
중앙빔과 그 좌우(가로), 상하(세로) 이웃빔의 3점을 이용해 각 축에서
독립적으로 포물선의 정점(peak) 위치를 보간한다.

한 축에서 등간격 $s$로 놓인 세 표본의 전력을
$y_L$(왼쪽/아래), $y_C$(중앙), $y_R$(오른쪽/위)라 하면, 이 세 점을 지나는
포물선의 정점은 중앙 표본 기준 정규화 오프셋

$$
\delta = \frac{1}{2}\,\frac{y_L - y_R}{\,y_L - 2y_C + y_R\,}\in[-1,1]
$$

로 주어진다(2차 보간 vertex 공식). 따라서

$$
\hat\theta_{az} = \phi_C + \delta_{az}\, s,
\qquad
\hat\theta_{el} = \psi_C + \delta_{el}\, s,
$$

- 가로축: $(y_L,y_C,y_R) = (A_W^2, A_C^2, A_E^2)$ → $\delta_{az}$
- 세로축: $(y_L,y_C,y_R) = (A_S^2, A_C^2, A_N^2)$ → $\delta_{el}$

코드의 `parabolic_offset()` 가 이 식을 구현하며, 분모가 0에 가깝거나
($\,y_L-2y_C+y_R\to0$, 즉 오목성 소실) 발산할 때를 대비해 $\delta$를
$[-1,1]$ 로 포화시킨다.

### 유도 (vertex 공식)
정점이 $\delta$인 포물선 $y(u)=a(u-\delta)^2+c$ 에 $u=-1,0,1$ 을 대입하면

$$
y_L=y(-1),\; y_C=y(0),\; y_R=y(1)
\;\Rightarrow\;
y_L-y_R = -4a\delta,\quad y_L-2y_C+y_R = 2a,
$$

두 식을 나누면 위의 $\delta$ 공식이 얻어진다.

### 특성
- **장점**: 빔당 단 3점으로 **서브-격자 연속 추정**, 가로/세로 분리로 계산
  단순, 단일 표적·고 SNR에서 매우 정밀(모노펄스의 진폭비교 추정과 동류).
- **모델오차(bias)**: 실제 빔 형상은 정확한 포물선이 아니라 sinc²(또는
  가우시안 근사)이므로, 표적이 빔 격자에서 멀수록(피크가 표본 구간 밖에
  가까울수록) 보간 편향이 커진다. 표적이 중앙빔 ±$s$ 이내일 때 가장 정확.
- **잡음 민감도**: 분모 $y_L-2y_C+y_R$ 가 작아지면(평탄/저SNR) 추정이
  불안정 → 포화 처리 필요. 전력 대신 **dB(로그) 값**으로 보간하면 가우시안
  빔에 더 잘 맞아 편향이 줄기도 한다(확장 옵션).
- **참고**: 가우시안 빔 가정 시 정확한 3점 추정식

$$
\delta = \frac{1}{2}\,
\frac{\ln y_L - \ln y_R}{\ln y_L - 2\ln y_C + \ln y_R}
$$

  (로그영역 포물선)이 흔히 사용된다. 본 구현은 선형전력 포물선을 채택.

### 용도
단일 표적의 정밀 각도추정(모노펄스 대체/보완), 빔 격자가 표적을 둘러싼 경우.

---

## 4. 비교 요약

| 항목 | Max-Amplitude | Centroiding | Parabolic Peak |
|------|---------------|-------------|----------------|
| 추정 해상도 | 격자 $s$ 로 양자화 | 연속(서브격자) | 연속(서브격자) |
| 사용 빔 수 | 1 (최대빔) | 전 빔(9) | 축당 3 (십자 5빔) |
| 계산량 | 최소 | 소 | 소 |
| 고 SNR 정밀도 | 낮음 | 중~높음 | 높음 |
| 저 SNR 강건성 | 높음 | 중 | 낮음(분모 불안정) |
| 주 오차원인 | 양자화 | 가장자리 편향 | 빔모델 불일치 편향 |
| 적합 상황 | 조대추정/탐지 | 중앙부·고SNR | 단일표적 정밀추정 |

### 이론 한계 (CRLB 관점)
각도추정 정확도의 하한은 Cramér–Rao Bound 로

$$
\sigma_\theta \ \gtrsim\ \frac{\theta_{3\text{dB}}}{k_m\sqrt{\text{SNR}}}
$$

형태로, 빔폭에 비례하고 $\sqrt{\text{SNR}}$ 에 반비례한다($k_m$: 빔
기울기 상수). Centroiding·Parabolic 은 고 SNR에서 이 한계에 접근하는 반면,
Max-Amplitude 는 격자 양자화로 SNR을 높여도 $s/\sqrt{12}$ 바닥에서 포화한다
(코드 Fig5: SNR 대비 RMSE 곡선에서 확인 가능).

---

## 5. 다중표적에서의 거동

세 추정기 모두 **단일 표적**을 전제로 한 단일 빔 클러스터 추정기다. 한
클러스터 안에 표적이 둘 이상 들어오면

- Max-Amplitude: 더 강한 표적의 빔으로 락(lock) → 약한 표적은 무시,
- Centroiding: 두 표적의 전력가중 중간점으로 끌림 → 양쪽 모두 편향,
- Parabolic: 합성된 단일 피크를 보간 → 표적 간격이 빔폭에 근접하면 붕괴.

따라서 표적 개수↑·표적 간격↓ 에 따라 RMSE가 급격히 증가한다(코드 Fig7,
Fig8). 다중표적 분해를 위해서는 클러스터 스캐닝, 반복 차감(CLEAN),
초분해능(MUSIC/ML) 등 별도 기법이 필요하다.

---

## 참고 문헌(개념)
- M. I. Skolnik, *Introduction to Radar Systems* — 안테나 패턴, 각도추정.
- S. M. Sherman, *Monopulse Principles and Techniques* — 진폭비교/모노펄스.
- D. K. Barton, *Radar System Analysis and Modeling* — 각도정확도/CRLB.
- 신호처리 일반: 3점 포물선 보간(quadratic peak interpolation).
