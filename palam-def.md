# 各パラメータの定義

|パラメータ名 | 定義 | 値域 | 用途 | 収集先 | 収集方法 |
| :--: | :-- | :-- | :-- | :-- | :-- |
 Estimated| プロジェクトが想定した開発工数(人日)| 1以上の整数 | 変数αの算出 | タスク管理ツール(構成管理ツール) | ソースコード開発に関する全タスクの工数を合計する. 収集できなかった場合, 構成管理ツールの変更保存ログのタイムスタンプから算出する. |
| Actual | 実際の開発工数(人日)| 1以上の実数| モデル構築手順proposed-step3のActual | 構成管理ツール | 構成管理ツールの変更保存ログのタイムスタンプから算出する. |
| FPTrial | FP試算法で算出した開発システムのFP値| 1以上の整数| モデル構築手順proposed-step2のSize | ソースコード | ソースコード中の関数からFP算出用のパラメータ集計と計算を行う |
| FPApproxi| FP概算法で算出した開発システムのFP値| 1以上の整数| モデル構築手順proposed-step2のSize(予備) | ソースコード | ソースコード中の関数からFP算出用のパラメータ集計と計算を行う |
| Requirement | ユーザからの全要求数 | 0以上の整数 | 工数変動要因 | 議事録 | 要件定義段階の議事録からユーザの要求数を集計する． |
| LikelySys | 開発システムに類似するシステムの調査件数 | 0以上の整数| 工数変動要因 | 議事録 | 開発予定のシステムに類似したシステムの調査件数を集計する． |
| Type| プロジェクトの開発種別| 新規開発はnew機能追加はadd再開発はre_develop| 工数変動要因  | 議事録 | 要件定義段階の議事録から開発システムの開発種別を確認する． |
| IsFormat| コードの記述例，フォーマットの有無| 0/1| 工数変動要因  | ソースコード | 複数ソースコードを比較し,表記の統一などが行われているかを確認する. |
| MeetSchedule| 作業開始から発表会，ヒアリングまでの日数| 1以上の整数| 工数変動要因 | タスク管理ツール | ソースコード開発予定日から成果発表の日程までの日数を算出する． |
| HearlingNum| ユーザへの総ヒアリング回数| 0以上の整数| 工数変動要因 | 議事録 | 要件定義段階の議事録からヒアリングの総実施回数を集計する． |
| AddiJobTime| 講義外の総会議時間(時)| 0以上の整数| 工数変動要因 | 議事録 | 設計段階までの全議事録から講義時間以外の会議時間を集計する． |
| Num| プロジェクト参加者数(人)| 1以上の整数| 工数変動要因 | 議事録 | プロジェクト初期段階の議事録から参加メンバー数を集計する． |
| StakeHolder| 想定ステークホルダー数(人)| 1以上の整数| 工数変動要因 | 議事録 | 要件定義段階の議事録からステークホルダー数を集計する． |
| WorkPasona| 1タスクあたりの同時参加メンバー数平均 | 1以上の実数 | 工数変動要因 | タスク管理ツール | 全タスクの担当者数を集計し，その平均値を有効数字3桁で算出する． |
| Importance| 各タスクの優先度を低い場合1，普通の場合2，高い場合3，緊急の場合4，今すぐ行わなければならない場合5としたとき1タスクの同時担当者数平均| 1以上の実数| 工数変動要因 | タスク管理ツール | 全タクスの重要度を設定し，¥¥その平均値を有効数字3桁で算出する． |
| IdentifyProb| 開発チーム内で発見した仕様上の問題点の洗い出し数| 0以上の整数| 工数変動要因 | 議事録 | 設計段階までの全議事録からプロジェクト中で洗い出した¥¥開発における問題点を集計する． | 
| NewTech| 開発メンバーの過半数が習得していない技術を利用しているか| 0/1| 工数変動要因 | 議事録 | 設計段階の議事録からメンバーの過半数が未習得の技術を¥¥利用するかどうか確認する． |
| Workspace | 作業スペースの確保度合い| 明らかに狭い場合1やや狭い場合2普通の広さがある場合3十分に広い場合4 | 工数変動要因 | 議事録 | 全議事録から会議実施場所情報を集計し，実際の状況を¥¥確認して数値を設定する． |
| SimilarProject| 類似プロジェクトの参加経験者有無| 0/1| 工数変動要因 | 議事録 | 要件定義段階の議事録から，プロジェクト参加者に¥¥類似プロジェクト参加経験があるか確認する． |
| PMTool| プロジェクト管理ツールの利用有無| 0/1| 工数変動要因 | 議事録 | 全議事録からプロジェクト管理ツールの利用有無を確認する． |
| CMTool| 構成管理ツールの利用有無| 0/1| 工数変動要因 | 議事録 | 全議事録から構成管理ツールの利用有無を確認する． | 
| DSTool| 設計支援ツールの利用有無| 0/1| 工数変動要因  | 議事録 | 全議事録から設計支援ツールの利用有無を確認する． |
| DCTool| ドキュメント作成ツールの利用有無| 0/1| 工数変動要因 | 議事録 | 全議事録からドキュメント作成ツールの利用有無を確認する． | 
| DTTool| デバッグ・テストツールの利用有無| 0/1| 工数変動要因 | 議事録 | 全議事録からデバッグ・テストツールの利用有無を確認する． | 
| CASETool| CASEツールの利用有無| 0/1| 工数変動要因 | 議事録 | 全議事録からCASEツールの利用有無を確認する． |
| CodeGenerator| コード自動生成ツールの利用有無| 0/1| 工数変動要因 | 議事録 | 全議事録からコードジェネレータの利用有無を確認する． | 
| DevMethodology| 開発方法論の利用有無| 0/1| 工数変動要因 | 議事録 | 全議事録から開発方法論の利用有無を確認する． |
| UserParticipation| ユーザの要求仕様関与度合い| 開発メンバーが全て作成する場合0ユーザが要求をすべて洗い出す場合1ユーザの要求仕様に修正を行う場合2すべてユーザが作成する場合3 | 工数変動要因 | 議事録 | 全議事録からユーザの¥¥関与度合いを確認する． |
| UserAgreement| 開発物に対するユーザ承認の有無| 0/1| 工数変動要因 | 議事録 | 全議事録からユーザの開発物に対する承認の有無を確認する． |
| MandatoryControl| 開発物に対する法的規制の有無| 0/1| 工数変動要因 | 議事録 | 全議事録から開発物に対する法的規制の有無を確認する． |
| Reliability| 要求される信頼性の度合い| 設定されていない場合0信頼性に関する項目が10件未満の場合110件以上の場合2 | 工数変動要因 | 議事録 | 設計工程の議事録からシステムの非機能要件に関する¥¥記述の量を集計する． |
| UC | ユースケース数| 0以上の整数| 工数変動要因 | 議事録 | 要件定義・設計工程の議事録から開発物のユースケース数を集計する． |
| Architecture| システムアーキテクチャの種類| アーキテクチャの名称| 工数変動要因 | 議事録 | 設計工程の議事録からシステムアーキテクチャの種類名を収集する． |
| Language| 主要使用言語| 言語名| 工数変動要因  | ソースコード | ソースコードの記述方法から利用する主要言語の確認を行う. |
| BizExp| 開発システムが利用される業務の従事経験の有無| 0/1| 工数変動要因 | 議事録 | プロジェクト初期の議事録から参加メンバーの開発対象領域の¥¥業務従事経験の有無を収集する． |
| DesignExp| システム設計経験の有無| 0/1| 工数変動要因 | 議事録 | プロジェクト初期の議事録から参加メンバーのシステム¥¥設計経験の有無を収集する． |
| LangExp| 使用言語の使用経験の有無| 0/1| 工数変動要因 | 議事録 | プロジェクト初期の議事録から参加メンバーの開発言語の¥¥使用経験の有無を収集する． |
| PlatformExp| 開発環境の使用経験の有無| 0/1| 工数変動要因 | 議事録 | プロジェクト初期の議事録から参加メンバーの開発¥¥プラットフォームの使用経験の有無を収集する． |
| PMExp| プロジェクトマネージャ経験者の有無| 0/1| 工数変動要因 | 議事録 | プロジェクト初期の議事録から参加メンバーのマネジメント¥¥経験の有無を収集する． |
| Grade| メンバーの学年平均(大学院生の場合学年+4とする)| 1以上の実数| 工数変動要因 | 議事録 | プロジェクト初期の議事録から参加メンバーの学年を¥¥¥ref{collected-metrix}の定義に従い計測し，平均値を算出する． |
| ProExp| システム開発経験者の参加率(歩合)| 0以上の実数| 工数変動要因 | 議事録 | プロジェクト初期の議事録から参加メンバーのシステム¥¥開発経験者の数を計測し，参加割合を算出する． |
| PBLExp| PBL経験者の参加率(歩合)| 0以上の実数| 工数変動要因 | 議事録 | プロジェクト初期の議事録から参加メンバーのシステム¥¥開発経験者の数を計測し，参加割合を算出する． |
| purpose| 開発物の開発目的が明確化されているか否か| 0/1| 工数変動要因 | 議事録 | 要件定義段階の議事録からプロジェクトの目的が記録¥¥されているか確認する． |  
| ConcuTask| 一人あたりの担当タスク数平均| 1以上の実数| 工数変動要因 | タスク管理ツール | 一人あたりの同時担当タスク数を集計し，その平均値を¥¥有効数字3桁で算出する． |
| ConcuProj| 一人あたりの同時参加プロジェクト数平均| 1以上の実数| 工数変動要因 | 議事録 | プロジェクト初期の議事録から参加メンバーの同時参加プロジェクト数を集計し，平均値を算出する． |