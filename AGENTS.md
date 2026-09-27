# FrontierPhysics agent 入口

固定共用治理： https://github.com/lijiabao1998/FrontierLab-Governance/tree/07d2b13051b83215182e411e1612f92f1912d8fb
讀其 AGENTS.md、RESEARCH_PROTOCOL.md、EVIDENCE_POLICY.md、SAFETY.md，再讀本庫 README、STATUS、VALIDATION 及題卡。

本次初始化後所有 agent 均走 `<agent>/PHYS-xxx-<topic>` 分支，作者不自合。owner 明確批准後才由 Integrator 合併；CODEOWNERS 不代表伺服器端分支保護已配置。fetch 後記 base SHA；一工作目錄一寫入者，先查領題與 PR。

每輪 start → 本輪四路聯網搜尋並讀原始來源 → 凍結驗收及預算 → admit → 重現 → 探索 → 獨立 verifier／skeptic → 留失敗 → PR。不能重用上一輪搜尋冒充本輪。無網路、缺資料授權、解答未核實或超預算就停。已確認外部同範圍解答立即標 COMPLETED_EXTERNAL，不冒領成果。

物理特別規則：明列假設、單位、近似、邊界、守恆和已知極限；分開模型誤差／數值誤差／統計誤差。精確極限、有限尺寸和實驗資料不能互相替代。未觀測現象只標預測，不宣稱實際發現。

研究檔案放 problems/<ID>/{experiments,proofs,results}/；runs/<round>/round.json 記參數、seed、環境、資料與結果 hash、執行命令、未做項。預設30分鐘／0美元／100次，有界結束。第一條主線 PHYS-001 基線；其他題先留佇列。
