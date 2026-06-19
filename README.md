# Radar Beam / Target Position Estimation Simulation

송신빔 1개 + 수신빔 9개(중앙 1 + 외곽 8, 3×3 사각배치) 구성의 레이더 빔
시뮬레이션과 표적 위치추정 알고리즘 성능분석 MATLAB 코드.

## 실행
```matlab
>> radar_beam_sim
```
MATLAB 기본 기능만 사용하며 별도 툴박스가 필요 없습니다
(`sinc` 는 로컬함수로 직접 구현).

## 모델 개요
- **도메인**: 각도공간(azimuth / elevation, deg)
- **빔패턴**: 직사각 개구면 sinc 패턴 — 안테나 크기(D)와 파장(λ)으로 빔폭 결정
  - one-way 전력패턴 `G = (sinc(D/λ·sin θ))²`, 송·수신 모두 적용
- **신호모델**: 빔별 수신전력 `P ∝ G_tx(표적) · G_rx,beam(표적) · RCS` (2-way),
  기준 SNR 기반 열잡음(복소 가우시안) 추가 후 포락선 진폭 측정
- **표적**: `[az, el, RCS]`, 다중표적은 독립 산란체(임의 위상)로 전압 합산

## 빔 배치 (조정 가능)
- 중앙 1 + 외곽 8 = 3×3 격자
- `cfg.beam_spacing` 으로 빔 간격(deg) 조정
- `make_rx_beams()` 의 `beam` 행렬을 직접 수정하면 임의 배치 가능
  (단, parabolic peak 는 중앙/상하좌우 십자 빔을 사용)

## 추정 알고리즘
1. **max-amp**: 최대 진폭 빔의 지향각
2. **centroiding**: 전 빔 전력가중 무게중심
3. **peak estimation**: 가로(W·C·E)/세로(S·C·N) 3점 포물선 보간으로 피크 위치 추정

## 분석/출력 (Figure)
- Fig1: 빔 배치, 송·수신 −3dB 빔패턴, 표적
- Fig2: 빔별 측정진폭 + 3종 추정결과
- Fig3: 단일표적 추정오차 산점도
- Fig4: 추정오차 크기 히스토그램
- Fig5: SNR 대비 RMSE
- Fig6: 보어사이트 오프셋 대비 RMSE
- Fig7: 다중표적 — 표적개수 대비 RMSE(최강표적 기준)
- Fig8: 다중표적 — 2표적 간격 대비 RMSE

## 주요 파라미터 (`cfg`)
| 항목 | 변수 | 기본값 |
|------|------|--------|
| 주파수 | `cfg.f` | 10 GHz |
| 송/수신 개구면 | `cfg.Dtx_az/el`, `cfg.Drx_az/el` | 0.30 m |
| 빔 간격 | `cfg.beam_spacing` | 1.0 deg |
| 기준 SNR | `cfg.SNR0_dB` | 30 dB |
| 몬테카를로 횟수 | `cfg.MC_N` | 3000 |

## 참고 / 한계
단일 빔클러스터(9빔)로는 한 클러스터 내 다중표적의 완전 분리가 원리적으로
제한됩니다. Fig7/Fig8 은 표적 개수·간격 증가에 따른 추정오차 열화를
정량적으로 보여줍니다.
