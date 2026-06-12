---
name: clinical-km-survival-analysis
description: Kaplan-Meier 생존분석 전용 스킬. 단일/다군 생존곡선, log-rank test, median survival(95% CI), number-at-risk table, 그리고 그룹 변수에 대한 unadjusted HR(univariate Cox)을 NEJM 스타일로 산출한다. 사용자가 Excel 파일을 업로드하고 "생존분석", "KM curve", "Kaplan-Meier", "생존곡선", "log-rank", "median survival"을 요청할 때 트리거. 다변량 보정(adjusted HR)·VIF·PH 가정 검정 같은 Cox 모델링이 필요하면 clinical-cox-regression 스킬을 사용하라.
---

# Kaplan-Meier Survival Analysis

비보정(unadjusted) 생존 비교 전용. KM 곡선 + log-rank + median survival + number-at-risk + 그룹 unadjusted HR.
다변량 보정·진단(VIF, Schoenfeld PH)은 `clinical-cox-regression` 스킬에서 수행한다.

## Workflow

### Step 1: 변수 정보 수집
사용자에게 확인:
1. **Time-to-event**: 추적기간 컬럼명 (days)
2. **Event**: 사건 발생 컬럼명 (1=event, 0=censored)
3. **Group**: 비교군 컬럼명 (선택)
4. **X축 간격**: 30일(단기) 또는 365일(장기)

### Step 2: KM 분석 실행
```bash
pip install lifelines openpyxl matplotlib pandas numpy --quiet --break-system-packages

python scripts/km_analysis.py --file <filepath> --time <time_col> --event <event_col> \
  --group <group_col> --interval <30|365> --output result \
  --labels 0:Control 1:Treatment
```
**Arguments**: `--file`, `--time`, `--event`, `--group`(선택), `--interval`(30 또는 365), `--output`, `--labels`(value:label 형식, 선택)

### Step 3: 결과 제시
1. KM 곡선 이미지 (number-at-risk table 포함)
2. 통계 요약표:
   - 군별 N, event 수
   - Median survival (95% CI)
   - Log-rank p-value
   - **Unadjusted HR (95% CI), p-value** — 그룹 변수에 대한 univariate Cox로 계산

## Output Format

| Group | N | Events | Median Survival (95% CI) |
|-------|---|--------|--------------------------|
| Control | 150 | 45 | 365 days (280-450) |
| Treatment | 148 | 32 | 520 days (410-NR) |

**Log-rank test**: χ² = 4.52, p = 0.033
**Unadjusted HR**: 0.68 (95% CI: 0.48-0.96), p = 0.028

## Notes
- Event coding: 1 = event, 0 = censored
- Time unit: days (입력값 그대로 표시)
- Graph style: NEJM colors, Arial, top/right spine 제거, number-at-risk table 필수
- **Unadjusted HR**은 단변량(그룹 1개) Cox이며, 교란보정이 아니다. 보정된 HR이 필요하면 clinical-cox-regression 사용.
