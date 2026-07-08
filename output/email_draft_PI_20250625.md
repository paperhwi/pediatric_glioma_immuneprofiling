# PI 보고 이메일 초안 (2025-06-25)

---

교수님 안녕하십니까,

소아교종 면역 생태형 원고 진행 상황을 보고드립니다.

## 1. 축소한 부분 (Cross-validation 구조 개편)

기존 7개였던 Cross-validation 섹션을 4개로 정리했습니다.

- **유지 (4개)**: Grabowski 모듈 검증(CV1), CD276 직교성 검증(CV2), Bagaev/PROGENy 시그니처(CV3), BRAF V600E LGG 재현 + TIS-PFS 다중비교 보정(CV4)
- **제거/흡수 (3개)**: 기존 CV5(Intermediate 생태형 생물학적 특성)는 Results 본문에 흡수, CV6(면역 시그니처 간 상관관계)은 CV1에 통합, CV7(생태형별 세포 비율 비교)은 이미 3.3절에서 충분히 기술되어 삭제

**축소 이유**: Reviewer가 "cross-validation이 과도하게 길다"는 지적을 할 가능성이 높고, 핵심 메시지(brain-tuned deconvolution → 3 ecotypes → 예후 연관)에 집중하는 것이 Cancers 투고 전략에 더 부합한다고 판단했습니다. 삭제된 내용은 모두 다른 섹션에서 이미 커버되거나, supplementary로 이동 가능합니다.

## 2. 6월 25일 추가 분석 내용

### (1) EOR(절제 범위) 민감도 분석
- **분석 내용**: 기존 Cox multivariable 모델에 EOR(GTR vs. non-GTR)을 공변량으로 추가
- **결과**: 생태형 HR이 EOR 보정 후에도 안정적 (Inflamed HR=0.50, p=0.027). EOR 자체는 유의하지 않음 (LR test p=0.58)
- **추가 이유**: Reviewer가 "수술 범위가 교란변수가 아닌가" 질문할 가능성이 높아 선제적으로 대응. Results 3.6, Limitations, Discussion에 반영 완료

### (2) CAR-T 타겟 발현 분석
- **분석 내용**: 소아교종 CAR-T 임상시험 주요 타겟 5개(B4GALNT1/GD2, ST8SIA1/GD3, IL13RA2, ERBB2/HER2, CD276/B7-H3)의 생태형별 발현량 비교 (Kruskal-Wallis + Dunn's post-hoc)
- **핵심 발견**: B4GALNT1(GD2 합성효소)이 생태형 간 유의한 차이 없음 (KW p=0.26) → **Immune-desert에서도 checkpoint blockade가 실패하는 환경에서 GD2-CAR T가 대안이 될 수 있음**을 시사
- **추가 이유**: Translational relevance 강화. 단순 면역 프로파일링을 넘어 치료 전략 제안까지 확장하여 논문의 임상적 가치를 높임

### (3) DESeq2 차별 발현 유전자 분석
- **분석 내용**: Inflamed vs. Immune-desert 간 DESeq2 (design: ~ cohort_group + ecotype, |log2FC| >= 1, padj < 0.05)
- **결과**: Inflamed-up 2,339개 (HLA-DRA +3.3, CXCL10 +4.3, IFNG +4.0 등 면역 시그니처), Desert-up 14,398개 (대부분 non-coding RNA)
- **핵심 발견**:
  - Checkpoint DEGs 확인: HAVCR2 (+1.8), PDCD1 (+2.1), IDO1 (+2.5) — LAG3/TIGIT는 역치 미달
  - B4GALNT1 non-DEG 재확인 (log2FC=-0.73) → CAR-T 분석과 일관
  - **VEGFA가 Desert에서 상향 발현** (log2FC=-1.19) → 항혈관신생 병용 전략의 근거
- **추가 이유**: ssGSEA 기반 생태형 분류가 실제 transcriptome-wide 차이를 반영하는지 독립적으로 검증. 또한 VEGFA 발견으로 "Desert 생태형에 대해서는 checkpoint 대신 항혈관신생 + CAR-T 조합"이라는 구체적 치료 전략 제안이 가능해짐

## 3. 향후 분석 방향 — 선택지

현재 원고는 Methods 15개 섹션, Results 11개 섹션, CV 4개로 구성되어 있으며, 이미 상당한 분량입니다. 다음 두 가지 방향 중 어느 쪽이 좋으실지 의견 여쭙겠습니다.

### Option A: 현재 결과로 마무리 (축소 방향)
- Supplementary figure 번호 정리 (S7, S_CART, S_volcano → 순차 번호)
- Figure legend 최종 점검
- Reference 정리 및 Cancers 양식 맞춤
- **장점**: 빠른 투고 가능, 이미 충분한 분석량
- **예상 소요**: 1-2일

### Option B: 분석 확장 후 투고 (확장 방향)
가능한 추가 분석:
1. **GSEA/Pathway enrichment**: DESeq2 결과에 대해 Gene Ontology, KEGG, Hallmark pathway enrichment 수행 → Desert-up 유전자의 생물학적 의미 규명
2. **Single-cell 검증**: Filbin et al. (2018) scRNA-seq 데이터로 우리 bulk 생태형 시그니처가 single-cell 수준에서도 재현되는지 확인
3. **생존 분석 세분화**: 생태형 × 분자아형(H3K27M, BRAF 등) 교호작용 분석
4. **면역세포 조성 비교**: CIBERSORTx absolute mode로 생태형별 세포 비율 정량 비교

- **장점**: 논문 depth 증가, higher-impact journal 가능성
- **단점**: 추가 2-4주 소요, Cancers 외 journal 재고려 필요할 수 있음

---

교수님 의견 부탁드리겠습니다.

감사합니다.
