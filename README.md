# PBL向け工数見積り手法の利用方法

## ※前提として

必要となるPBL実績データを収集してください．
工数変動要因は順不同でも構いませんが，Estimated，Actual，FPTrial，FPApproxiは必ずこの順でデータの先頭に配置してください．

収集した過去のPBL実績データは通常のPBLと短期集中演習型のPBLとに分割し，それぞれ”normal”，”mini”という列名でデータフレームに記録してください．
見積り対象PBLのデータはこの分割は不要ですが，上記の順番でデータを記録してください．

## 特定のPBLを対象に見積りを行う場合

`CalcPBLManHours(TeachingData, TestData, PBLType = "normal")`
を利用する．

`TeachingData`は過去のPBL実績データ．
`TestData`は見積り対象のPBLデータ．
`PBLType`は見積り対象PBLの演習形式が通常演習か短期集中演習かを指定する．通常演習の場合`normal`，短期集中演習の場合`mini`を代入する．無記入の場合`normal`として動作する．

## 過去の実績データのみで見積り実験を行う場合

`CalcLOOCV(Dataset)`を利用する．

`Dataset`は過去の実績データ．ただし代入する際にはnormal”，”mini”を含めたデータフレームではなく，`Dataset$normal`のように対象を指定すること．

## 見積り結果の相対誤差を計算する場合

`CalcRelativeError(experiment, calculated)`を利用する．
計算結果は%形式で出力される

`experiment `は見積り結果の値(または値の配列)．
`calculated `は見積り対象の実績値(または実績値の配列)．

## 二乗平均平方根を計算する場合

`RMS(data)`を利用する．

`data`は計算対象の数値の配列．

## 相対誤差の分布を図に出力する場合

`plotError(error)`を利用する．
図中の日本語はヒラギノ角ゴProフォントを利用している．

`error`は見積り結果の相対誤差の配列．