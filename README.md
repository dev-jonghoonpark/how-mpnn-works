# How MPNN Works

**MPNN — 원자들이 이웃에게 말을 건네는 신경망.** Gilmer 외(Google Brain · DeepMind, ICML 2017)의 「Neural Message Passing for Quantum Chemistry」를 브라우저에서 직접 실험하며 이해하는 인터랙티브 교육 자료입니다. 따로 나온 그래프 신경망 여덟 개가 메시지 M · 갱신 U · 판독 R 세 함수의 선택만 다르다는 통일, 그리고 그 틀 안에서 찾은 최고 모델(엣지 네트워크 + set2set)이 QM9 13개 성질 모두에서 최고 성능, 11개에서 화학적 정확도에 도달한 과정을 따라갑니다.

🔗 **[데모 바로가기](https://dev-jonghoonpark.github.io/how-mpnn-works/)**

## 다루는 내용

1. **병목은 정확도가 아니라 비용** · DFT는 QM9 분자 하나에 코어 1시간, 신경망 추론은 30만 배 빠르다 — 후보 분자 수·코어 수로 보는 스크리닝 시간 계산기
2. **분자를 그래프로 적기** · 논문 표 1의 원자 특징, 엣지 표현 세 가지(화학 그래프 / 거리 구간 14개 알파벳 / 실수 거리 5차원)를 분자 위에서 확인
3. **식 (1)~(3)** · 장난감 MPNN(d = 4)을 실제로 돌려 원자를 고르고 t를 올리며 메시지 합과 갱신 값을 한 단계씩 보기, 정보가 그래프 거리 t 안으로 퍼지는 모습
4. **순열 불변** · 원자 번호를 섞으면 합 판독은 그대로(차이 ~10⁻¹⁷), 이어 붙이기 판독은 바뀜
5. **여덟 개의 모델, 하나의 틀** · Duvenaud · GG-NN · 상호작용 네트워크 · Kearnes · DTNN · 스펙트럼 · GCN · enn-s2s의 M, U, R
6. **메시지 함수와 엣지 표현** · 거리 구간 행렬(계단) vs 엣지 네트워크 A(e)(연속), 표 9 쌍 메시지 vs 엣지 네트워크, 표 10 거리 정보의 효과(2.57 → 0.98 → 0.68)
7. **가상 엣지와 마스터 노드** · T에 따른 도달 행렬, 모두 닿는 데 필요한 T와 메시지 수, 표 3
8. **판독: 합 vs set2set** · 합이 같은 세 집합을 어텐션 판독으로 구별하기
9. **타워** · n²d² → n²d²/k 비용 계산기, 실측 약 2배, 표 4
10. **QM9 결과** · 표 2를 성질별 로그 눈금 점 그래프로(오차 비율 1 = 화학적 정확도), 표 6 학습 데이터 크기에 따른 오차 비율
11. **브라우저에서 학습하는 MPNN** · 거리 3인 헤테로원자 쌍 세기를 T = 1, 2, 3으로 동시 학습 — T = 1은 기준선에 멈추고 T = 2부터 맞히는 이유(판독이 ⌈K/2⌉를 만든다)
12. **MPNN 이후** · 큰 그래프로의 일반화, 메시지 어텐션, SchNet · DimeNet · 그래프 트랜스포머
13. **퀴즈** 6문항

## 참고

- Gilmer, Schoenholz, Riley, Vinyals & Dahl (2017), *Neural Message Passing for Quantum Chemistry* — [arXiv:1704.01212](https://arxiv.org/abs/1704.01212)
- Li 외 (2016), *Gated Graph Sequence Neural Networks* — [arXiv:1511.05493](https://arxiv.org/abs/1511.05493)
- Vinyals, Bengio & Kudlur (2016), *Order Matters: Sequence to sequence for sets* — [arXiv:1511.06391](https://arxiv.org/abs/1511.06391)
- Faber 외 (2017), *Prediction errors of molecular machine learning models lower than hybrid DFT error* — [arXiv:1702.05532](https://arxiv.org/abs/1702.05532)

[how-ai-works](https://github.com/dev-jonghoonpark/how-ai-works) 시리즈의 일부입니다.
