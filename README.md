# FrontierPhysics
前沿物理

## 主線：可重現基線 → 物理一致性 → 可區分預測
第一輪先處理 **PHYS-001 的基線與驗證器設計**：凍結一個可負擔模型，重現能量平衡及小尺度結果，不直接宣稱解決三維湍流。其餘問題先排隊，避免十條線同時消耗算力。

初始題目是從近期研究整理的「尚未完成的一般問題／可研究缺口」。`OPEN` 只表示本次有界篩查未確認同範圍解答；不是全網不存在答案的證明。每輪仍須重新聯網。

| ID | 問題 | 優先級 | 第一個可做的小任務 |
|---|---|---|---|
| [PHYS-001](problems/PHYS-001/problem.json) | 三維湍流間歇性與異常尺度律 | A | 能量平衡、合成資料及結構函數基線 |
| [PHYS-002](problems/PHYS-002/problem.json) | 摻雜二維 Hubbard 模型的熱力學相圖 | B | 小晶格精確對角化與已知極限 |
| [PHYS-003](problems/PHYS-003/problem.json) | 受限多體模型的符號問題與無偏取樣效率 | A | 固定小 Hamiltonian，比較方差及成本 |
| [PHYS-004](problems/PHYS-004/problem.json) | 多體局域化的熱力學穩定性 | B | 有限尺寸、時間及無序敏感性 |
| [PHYS-005](problems/PHYS-005/problem.json) | 有限維理想玻璃轉變與鬆弛機制 | C | 區分模型預測與可觀测交叉轉變 |
| [PHYS-006](problems/PHYS-006/problem.json) | 中子星高密度物態方程的可識別性 | B | TOV 合成資料與先驗敏感性 |
| [PHYS-007](problems/PHYS-007/problem.json) | Hubble tension 的系統誤差與模型解釋 | B | 公開似然及誤差模型重現 |
| [PHYS-008](problems/PHYS-008/problem.json) | 暗物質身分及直接探測的訊號識別 | C | 背景模型、統計覆蓋與新結果核查 |
| [PHYS-009](problems/PHYS-009/problem.json) | 黑洞資訊編碼：受控模型到物理時空 | C | 重現一個模型並列出關鍵假設 |
| [PHYS-010](problems/PHYS-010/problem.json) | 暗能量是否隨時間演化 | B | 凍結資料組合、先驗與模型比較 |

A/B/C 是啟動成本及排序，不是解題概率。每張卡都有 statement、known_result、open_gap、第一步、evaluator、限制、完成條件及來源。物理題不因 benchmark 改善就整題完成。

## 每輪開始
讀 [AGENTS.md](AGENTS.md)、[STATUS.md](STATUS.md)、[VALIDATION.md](VALIDATION.md) 和固定治理版本。

```bash
# 治理庫與本庫放同一父目錄，治理 checkout GOVERNANCE.lock.json 的 commit
python3 ../FrontierLab-Governance/tools/frontier.py validate .
python3 ../FrontierLab-Governance/tools/frontier.py start . PHYS-001 --agent gpt
# 真正做本輪 general / discipline / solution / criticism 檢索，讀原文並填 round.json
python3 ../FrontierLab-Governance/tools/frontier.py admit . runs/<round-id>/round.json
```

同範圍問題已有獨立核實解答：立即標 `COMPLETED_EXTERNAL`、附作者及來源、停止重複探索。新模型或預印本自稱解答但有爭議：`CLAIMED_RESOLVED`，先查證。只改善某參數區域：保留父問題，更新局部結果。

## 驗證邊界
量綱、守恆、已知極限、數值收斂與資料相容性分別檢查，並記適用條件。toy model 不是自然界；有限尺寸不是熱力學極限；擬合成功不是唯一機制。每個子任務在動手前寫定參數、邊界、對照、誤差與停止條件。

## 目前狀態
本次建立題卡及流程，沒有執行原創物理實驗，也没有配置 GPU、常駐 agent 或付費 API。共用 CI 只驗紀錄完整性與每輪前置檢索紀錄；不能替代真正閱讀文獻和物理驗證。
