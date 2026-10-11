# PHYS-001 小 N 獨立元件驗證報告

本輪 **SYNTHETIC_COMPONENT_AUDIT_PASS**，42 個具名檢查全部通過、退出碼 **0**。這只驗證固定有限周期合成訊號的 DFT/FFT 與能量控制；原 r1 的整輪 **BASELINE_FAIL_DO_NOT_PROCEED_TO_NOVELTY** 保留，沒有真湍流發現宣稱。

核對來源是 [PR1](https://github.com/lijiabao1998/FrontierPhysics/pull/1) head `3740370457e32f64b22ac1bce8a8ea3467f198f8`，修正對照是 [PR2](https://github.com/lijiabao1998/FrontierPhysics/pull/2) head `51107b50e4c1a855c13ec1eb87381f480713f547`；fetch 後 main 為 `bca6b6c16efc5d2e04ff95cdfda5bea522160371`。治理鎖定 `07d2b13051b83215182e411e1612f92f1912d8fb`。父代理已登記本元件審計，與其他 PR 的結構函數/視窗修復範圍分開；沒有新的 GitHub 寫入。

計畫在執行前經另一個唯讀 session 審查、完成四路新查與本地 admit。協議 SHA256 為 `92928d4943e4b6cf59c2b63a3279ecdf528af2f42342fbc1fa15b52c7e6e8c81`；先行本地協議提交 `897adeb`。准入後未調整任何判準。題卡的物理問題比本輪元件範圍更廣。

固定 N=32/64/128，各自測 zero、DC、cosine、sine、Nyquist、mixed，另外有無 DC/Nyquist 的能量訊號。完全確定，沒有隨機 seed 或外部資料。內部模式是 A*cos(2π*k*j/N+phase)，mixed 的 (k,A,phase) 是 (3,1.25,0.3)、(5,0.75,-0.7)、(7,0.4,1.1)，DC=0.25、Nyquist=0.5；能量控制保留相同模式但 DC=Nyquist=0。單位均為無因次，周期為 1，沒有 PDE 或時間演化。

獨立實作直接求和 `X[k]=Σu[j]*(cos(2πjk/N)-i*sin(2πjk/N))`，使用 math.fsum，沒有借用 FFT。原碼只擷取經逐函式檢查的 fft、energy_budget AST，將 N 綁定到當前 fixture；原 module top-level、make_field、main 沒有執行。符號及 forward 無正規化約定由 [FFTW 定義](https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html) 與 [NumPy 公式](https://numpy.org/doc/stable/reference/routines.fft.html) 核對；未使用兩個套件的執行庫。

凍結的複係數與能量绝對誤差上限是 1e-10。正頻內部模式係數是 N*A*exp(i*phase)/2，sine 的正頻係數是 -i*N*A/2。平均 u² 的 Parseval 值是 Σ|X|²/N²。

| 實際指標 | 全部小 N / fixture 最大誤差 | 凍結上限 |
|---|---:|---:|
| DFT/FFT/解析複係數 | 1.3116551764792642e-12 | 1e-10 |
| Parseval 平均能量 | 4.218847493575595e-15 | 1e-10 |
| 解析平均能量 | 2.220446049250313e-16 | 1e-10 |

FFT r1/r2 均與直接 DFT 及解析 signed coefficients 一致。原相對能量容差 0.10 保持不變：每個 N 的乾淨輸入 PASS，1.5x 輸入在原參照振幅不變時 FAIL，膨脹殘差接近解析值 1.25。

符號共軛 mutation 被複係數檢查拒絕，即使 Parseval 能量仍相同；forward 除 N、Parseval 誤除 N 或 N³ 都被拒絕。將膨脹輸入同時配上 1.5x 參照振幅，能量 guard 會誤接受，因此負控檢查將此 mutation 拒絕。錯誤的 [1.0] 參照使乾淨訊號不通過，也被拒絕；不能只看到膨脹訊號被拒就稱 guard 有辨識力。缺失、失敗或空總檢查也都拒絕成功。所有 42 個具名檢查及 21 組完整 samples 保存。

既有基線與本輪數值分開：

- r1 儲存的 E3=1.5400285073651272 < 1.7、check=false、整輪 FAIL 完全保留，未重算。
- r2 amplitude guard 的乾淨殘差 4.090885079679787e-14、膨脹殘差 1.250000000000092 是 **INHERITED_NOT_RERUN**。修正使用真振幅已存在於 PR2；本輪的小 N 核對不是新修復或新成果。
- r1 既有 results/r1/hashes.txt 為結果 JSON 記錄 88dd5b...，實際固定 Git blob 是 bd9e814...。此歷史 checksum 失配在 INHERITED_INTEGRITY.json 中保留，沒有更寫舊 manifest。新元件 PASS 不代表原整輪或所有歷史證據已通過。

實際啟動 2026-10-11T18:23:51.510445Z，結束 18:23:51.845470Z。单 CPU affinity mask=1；子程序總 wall 0.172 秒，harness wall 0.047 秒、CPU 0.031 秒，$0。數學子程序在 18:23:04 已結束，兩計算不重疊。stderr 保留 Python 啟動路徑警告，原因未判定；完整新结果确已產出。

復跑在庫根目錄執行，使用新 tag：

```powershell
python -X utf8 problems/PHYS-001/experiments/independent_20261011/run_bounded.py --round runs/20261011T180848438131Z-gpt-PHYS-001/round.json --tag replay-01 --replay
```

預期退出碼 0。run_bounded.py 核對凍結 input hashes、限制單 CPU、最多30分鐘，並拒絕覆蓋既有 attempt。原始 launch、execution、stdout、stderr 與 result.json 位於 attempt-01/，完整實數輸入與複係數向量在 attempt-01/data/samples.json。

獨立唯讀審查檢查了協議、程式、來源 hash 及原始資料；並未第二次執行 DFT 或完整物理基線，不能宣稱第二次獨立實驗。自寫直接 DFT 與原 radix-2 FFT 是不同數值算法；作者已讀原碼，因此不是盲測。最終審查在 INDEPENDENT_REVIEW.md。

Git 初始文字轉換造成 3 個來源快照的提交 blob 與實際執行 bytes 不同，失配在 BYTE_PRESERVATION.json 保留。後續僅本輪子目錄使用 `* -text`、重新 stage 實際 bytes，最終 Git-blob hash 核對為交付要求。研究程式、協議、門檻與原始數值均未因此修改或重跑。

未跑 N4096 整輪、r2 校準、原 amplitude_guard.main、結構函數/視窗分析、DNS、Navier–Stokes、統計獨立性或收斂研究。模型誤差仍是合成訊號無法代表真湍流；數值誤差是表中有限浮點誤差，沒有統計估計。本輪 [Mukherjee 等來源](https://arxiv.org/abs/2607.26896) 只閱讀摘要與版本史，未作全文物理理論審計。PHYS-001 保持 OPEN。只本地提交，由父代理協調回庫，不推送、不合併、不部署。
