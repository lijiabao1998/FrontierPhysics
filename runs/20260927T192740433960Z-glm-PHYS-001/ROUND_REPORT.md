# ROUND REPORT｜PHYS-001 r2｜20260927T192740433960Z-glm-PHYS-001

- branch：`glm/PHYS-001-calib-r2`；pin `07d2b130`；verdict：**NO_RESOLUTION_FOUND**（primary 再核查：v1、無 erratum、理論缺口原文確認）
- 本輪判定：**循環 1 FAIL（如實記錄）→ 循環 2 CALIBRATED_BASELINE_PASS**——兩段式預登記的完整演示

## 方法（修正 r1 之門檻校準失誤）
兩段式：Stage A 校準（只用 Stage A 資料凍結帶）→ Stage B 確認（全新 seeds）。r1 的 E3 FAIL 未回頭改——以新流程取代。

## 循環 1（seeds A:1-8 / B:11-20）——FAIL，根因有價值
帶規則「均值±4σ」在估計器近確定性（跨 seed σ≈1e-6）下退化為**零寬度帶**（例：[0.6384, 0.6384]），Stage B 任何浮點級差異即出帶 → B1 FAIL。附帶量測：N 敏感性（2048→4096 斜率位移 0.07–0.18，β=5/3 窗 W1 為 0.107）；窗敏感性（三窗 spread：β=5/3 0.084、β=2 0.110、β=3 0.179——長窗與邊際 β 更敏感）。
**教訓：σ→0 時 ±kσ 帶退化；預登記規則必須含絕對下限。**

## 循環 2（seeds A:21-28 / B:31-40，帶=均值±max(4σ, 0.01)）——PASS
- Stage A 後凍結 bands2.json（例：β=5/3|W1 SF 帶 [0.6284, 0.6484]，中心 0.638≈2/3 ✓）
- Stage B：3β × 10 seeds × 3 窗 = **90 項量測 0 違規**；頻譜斜率全在帶內；能量殘差全帶內
- B4 分離度：三 β 之帶互不重疊（0.64 / ~1.0 / ~1.54）✓
- verdict：**CALIBRATED_BASELINE_PASS**（凍結帶於 Stage B 前寫入 bands2.json，時間戳可稽）

## 交付
`experiments/r2_calibration/{calibrate_and_confirm.py, bands.json(循環1), bands2.json(循環2)}`、`results/r2/{r2_calibration_results.json, r2_cycle2_results.json, hashes.txt, environment.txt}`、LITERATURE_MAP、本報告。

## 沒做成／限制
- Sabra/GOY shell model 未實作（時間盒）——guards 已規格化：能量收支 ⟨f·u⟩=⟨νk²|u|²⟩（容差預登記）、stationarity（窗口均值漂移檢定）、慣性區 k^-5/3 斜率帶（沿用本輪校準法）。
- 合成場仍非 NS 解；估計器僅二階。

## 下一輪最小下一步
1. 以循環 2 凍結帶為固定驗收，實作 Sabra（stdlib RK4，N殼=22，f 於 n=4-6）——三 guards 全綠才收數據。
2. 把帶規則（±max(4σ, floor)）寫入未來所有兩段式流程的模板。

## Remediation v2（Convergence Wave，Codex 二審回應）
- **P1 cycle-2 executable**：`calibrate_and_confirm_cycle2.py` 已入庫並實際重跑——重現 bands2.json＋r2_cycle2_results.json（Stage A seeds 21-28、Stage B 31-40、floor 規則、0/90 違規）。
- **P1 語義降級（Codex P1，採方案 B）**：整週期整數 Fourier 模態場之 seed 只變 phase，能量/頻譜/全域 S2 對 phase 不敏感（交叉項相消）——Stage A/B 是**確定性數值重現檢查**，不是統計確認；真統計泛化需變幅/譜擾動/獨立有限子窗（延後）。
- **P1 hashes**：verify_manifest.py，最後生成。
- **P2 B4**：三個凍結窗全部檢查分離（W1/W2/W3 全 true）。
- **P2 amplitude guard**：改用 make_field 真正返回之 amplitudes——unmodified 4.1e-14 PASS、1.5× 1.25 FAIL（results/r2/amplitude_guard_check.json）。
- **dsk verifier 歸因**：β=3 之 1.54 是**窗依賴統計量**（marginal log 修正），非 clean exponent、非 universal exponent——依 dsk 獨立審計（PR #3）記錄；其「finite-window saturation」結論僅部分吸收（β=3 log 修正保留），不宣稱 dsk 已排除 implementation error（其自身 review 未收斂）。
