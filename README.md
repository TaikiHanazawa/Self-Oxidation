# Approach

---

## 1. ゴール

- Sakuraba et al. (2021) の V-shaped（C-N-H）と，Shi et al. (2022) の δ¹⁵N（atm ~ 0‰, mantle ~ −5‰）の同時再現
- **crust は対象外**：crust の δ¹⁵N（~+5‰）は生物的 N 固定・堆積物蓄積・プレートテクトニクスによる後生的プロセスで形成されるため，集積モデルで再現すべき対象ではない

  → **atm（0‰）と mantle（−5‰）の 2 層を再現目標とする**
- 論文ベースで実装．物理的に妥当な説明が難しいパラメータは使わない．

---

## 2. Impactor シナリオ（単位は M_Earth）

### 2-1. 成長フェーズの区切りと GI 発数

**フェーズ一覧：**

| 質量範囲 | フェーズ |
|---|---|
| 0.01 → 0.10 | 暴走成長・寡占成長フェーズ（連続集積） |
| 0.10 → 0.90 | 巨大衝突フェーズ GI1–GI8（各 0.1 M_Earth，EC-like × 7 発 + CC event × 1 発） |
| 0.90 → 0.995 | GI9 = Theia（月形成 GI，0.095 M_Earth，EC-like） |
| 0.995 → 1.0 | late veneer（EC-like） |

**発数の算出：**
- 0.10 + 8 × 0.10 = 0.90（GI1–GI8）
- 0.90 + 0.095 = 0.995（GI9 = Theia，端数）
- → 合計 9 発，端数（0.095）は最後に配置

**初期条件の根拠**
- 初期揮発性元素は EC bulk 組成として MO 中に一様に設定し，コア未分化として扱う
- 1 au では oligarchic growth により `~10^26 g ≈ 0.017 M_Earth` の protoplanet が `~5×10^5 yr` で形成されるため，`0.01 M_Earth` を embryo-scale start point とみなすのは動力学的にも自然である（Kokubo & Ida 2000, *Icarus* 143, 15-27）
- 0.10 M_Earth をフェーズ境界とするのは，1 au 付近の terrestrial region で oligarchic growth により形成される protoplanet / planetary embryo の isolation mass が `~0.1 M_Earth` であり，その先は embryo 同士の giant impact stage に移るため（Kokubo & Ida 2012, *Prog. Theor. Exp. Phys.* 2012, 01A308）

**GI8 を CC event とする根拠**（0.80 → 0.90 M_Earth，地球質量の 80–90% 時点）：

- N 体計算では，CC（外部太陽系）物質の地球への到達は地球質量の 60–90% 形成後に集中
  - O'Brien et al. (2006): last 20–40%（Earth 60–80% 以降）
  - Rubie et al. (2015): 60–80% 以降で CC 的物質が到達
- Mo/Ru 同位体制約では CC 物質の混入は最後の 10–20%（Earth 80–90%）に特定
  - Kleine et al. (2020, Space Sci. Rev.); Budde et al. (2019, Nature Astron.)
- Izidoro et al. (2022) の N 体計算では Earth ~90% 時点に CC embryo 1 発が衝突
- 以上から地球質量 80–90% 時点（= GI8）が CC 物質到達の最も論文的根拠の強い位置

**GI9 = Theia（月形成 GI）の根拠：**

- Hf-W 年代から Theia 衝突は CAI 後 ~40–100 Myr（最良推定 ~50–62 Myr）
  - Kleine et al. (2009); Touboul et al. (2007)
- 衝突時の地球質量 ~90% M_Earth（= GI9 開始時点 0.90 M_Earth と整合）
- 質量 0.095 M_Earth は 0.1 M_Earth とほぼ等しく，Theia ~0.1 M_Earth という標準的想定と整合
  - Canup & Asphaug (2001, Nature)

---

### 2-2. Impactor 組成と δ¹⁵N

**δ¹⁵N の設定値：**

| 種別 | δ¹⁵N | 根拠 |
|---|---|---|
| EC-like | −30‰ | EC の δ¹⁵N 範囲は −45‰〜−15‰，代表値として採用（Shi et al. 2022 Fig.1; Grady et al. 1986） |
| CI-like | +40‰ | CI および CM コンドライトの実測 δ¹⁵N 範囲は +42〜+175‰（Shi et al. 2022 本文）；Shi et al. (2022) が accretion モデルでの入力値として +40‰ を採用（≈実測範囲の下限）；本モデルもこの値に従う |

#### 暴走・寡占成長フェーズ（0.01 → 0.10）：生 EC 組成

以下の数値はすべて Sakuraba et al. (2021) Table 1 の "EC model" 欄の値：

| 元素 | 濃度 |
|---|---|
| C | 4,000 ppm |
| N | 250 ppm |
| H | 400 ppm |
| δ¹⁵N | −30‰ |

**根拠：**
- 水スノーラインの圧力バンプが NC と CC 物質を空間的に分離し，内側リング（~1 au）は NC（非炭素質）材料のみから形成される（Izidoro et al. 2022）
- 多元素同位体（O, Ti, Cr, Ni, Mo, Ru）の質量収支から，地球の主成分は ~91% enstatite chondrite（EC）起源（Dauphas 2017, Nature）
- 上記 2 つの根拠から EC 組成を採用し，具体的数値に Sakuraba et al. (2021) Table 1 を使用

#### 巨大衝突フェーズ GI1–GI7（0.10 → 0.80）：分化済み EC 体の組成

| 元素 | 濃度 |
|---|---|
| C, N, H | 0.10 M_Earth 時点のモデル出力値を使用 |
| δ¹⁵N | −30‰（EC-like） |

**根拠：**
- 0.10 M_Earth の impactor は，それ自身が EC 素材から成長・コア分化した分化天体であり，生の EC 組成をそのまま適用することは不適切
- Sakuraba et al. (2021) が同様の方針を採用：*"we calculated the abundances of C, N, and H in protoplanets by running our model for the growth by planetesimal accretion in advance from 0.05 to 0.1 Earth masses"*
- 本モデルでは 0.01 → 0.10 の計算出力が自然にこの組成を与えるため自己整合的

#### 巨大衝突フェーズ GI8（0.80 → 0.90）：CC event（CI-like）

| 元素 | 濃度 | 出典 |
|---|---|---|
| C | 35,000 ppm | Sakuraba et al. (2021) Table 1 "CI model" |
| N | 1,500 ppm | Sakuraba et al. (2021) Table 1 "CI model"（Shi et al. のモデルでは ~1000 ppm；本モデルは Sakuraba 値を採用） |
| H | 6,900 ppm | Sakuraba et al. (2021) Table 1 "CI model" |
| δ¹⁵N | +40‰ | Shi et al. (2022) |

**根拠（CC event であることの正当化）：**
- GI8 が地球質量 80–90% 時点であることは Mo/Ru 同位体・N 体計算から支持される（Section 2-1 参照）
- 「δ¹⁵N = −5‰ の mantle」を作るには，EC（−30‰）だけでは不十分で，高 δ¹⁵N 物質を MO が液体の状態で受け取る必要がある
- Shi et al. (2022) が明示：CI 的物質は MO 中に投入され，後続の giant impactor の金属が MO 内の N（EC + CI 混合）をコアに分配する → 残留マントルに ¹⁵N が濃縮
- CC 的物質の GI 期混入は Grand Tack（Walsh et al. 2011）による動力学的散乱で説明可能
- 本モデルでは GI8 impactor 自体は `~0.1 M_Earth` 級の**分化済み CC embryo**として扱う．揮発性元素の bulk inventory は CI-like 値を用いる一方，`D_N`, `D_H`, `D_C` に入る金属相 alloy 組成には EC-like proxy を用いる．具体的には `x_S`, `x_Ni`, `x_Si` は EC representative values を保ち，`x_C` は 0.01→0.10 `M_Earth` precompute の出力コア組成から得た値を用いる．この proxy は分配係数計算にのみ用い，pre-existing な分化コアの揮発性濃度 `C_imp_met,i` や `δ¹⁵N_imp_met` を GI8 に別途与えることはしない．また self-oxidation に効く GI8 silicate redox は CC 固有値を新たに仮定せず，衝突時点の target mantle redox を用いる．

#### GI9 = Theia（0.90 → 0.995）：EC-like（分化済み）

| 元素 | 濃度 |
|---|---|
| C, N, H | GI1–GI7 と同じ（0.10 M_Earth 時点のモデル出力と同組成） |
| δ¹⁵N | −30‰（EC-like） |

**根拠：**
- Theia（月形成 impactor）は EC を主成分とする分化天体と考えられる
- 月の同位体組成（O, Ti, W）は地球との高い類似性を示し，Theia も地球と同じ内側リング（NC）起源であることを示唆（Pahlevan & Stevenson 2007; Dauphas 2017）

#### late veneer（0.995 → 1.0）：EC-like 組成

| 元素 | 濃度 | 出典 |
|---|---|---|
| C | 4,000 ppm | Sakuraba et al. (2021) Table 1 "EC model" と同値 |
| N | 250 ppm | 同上 |
| H | 400 ppm | 同上 |
| δ¹⁵N | −30‰ | — |

**根拠：**
- Ru 同位体が CC（外部太陽系）起源の late veneer を否定（Shi et al. 2022）
- Kuwahara et al. (2023) が EC の金属 Fe 降下によるマントル還元効果を明示
- late veneer の主効果は N 侵食（erosion 主体）であり，δ¹⁵N への影響は小さい．ただし EC-like N（−30‰）の投入が大気 δ¹⁵N を若干負に引く効果は計算で確認する（Sakuraba et al. 2021 の大気侵食方程式を使用）


---

## 3. Self-Oxidation

### 3-1. 概要

各衝突ステップで以下のチェーンを計算する：

```
impactor 質量・速度
  → melt pool 深さ
  → melt pool 底圧力 P_eq
  → self-oxidation と activity 比から更新された ΔIW_eq,deep（金属-ケイ酸塩平衡条件，Section 3-6A）→ D_N, D_H, D_C へ渡す
  → Fe³⁺/ΣFe(P_eq, T, ΔIW_eq)（Section 3-5）→ マントル Fe³⁺ 収支・表面 ΔIW（Section 3-6B）
```

**self-oxidation の圧力依存性**
- Zhang/Hirschmann の式は `P ≈ 0` から連続的に `Fe³⁺/ΣFe(P,T,fO₂)` を与えるため，本モデルでは明示的な on/off 閾値は設けない
- ただし Armstrong et al. (2019) の peridotitic melt 実験では，pressure-induced な Fe³⁺ 安定化が顕著になるのは `P > 10 GPa` からであり，これは self-oxidation が本格化する代表スケールの目安とみなす

---

### 3-2. Melt Volume の計算（de Vries et al. 2016）

**衝突角度補正を含む体積式**（Abramov et al. 2012; de Vries 2016 Eq.2）：
```
V_melt = (π/6) × k × E_m^(−3μ/2) × (ρ_p/ρ_t) × D_p³ × v_imp^(3μ) × sin^(2p_sin)(θ)
```

**スケーリング定数**（de Vries 2016 Table 1, Abramov et al. 2012 より）：

| 定数 | 値 | 備考 |
|---|---|---|
| k | 0.42 | — |
| μ | 0.56 | — |
| p_sin | ≈ 0.66 | sin(θ) にかかる指数．式 3μ/(2+μ) = 3×0.56/(2+0.56) = 0.656 を 0.66 に丸め（de Vries 2016 Table 1; de Vries 原文では γ と表記，ここでは Nakajima 2021 の質量比 γ との混同を避け p_sin と呼ぶ） |
| θ | 45°（中央値） | Nakajima et al. (2021): N 体計算（Rubie et al. 2015）の統計で 30°≤θ≤60° が典型的；Canup (2004): 月形成 GI も θ ≈ 45° |

**衝突速度**（de Vries 2016; Nakajima et al. 2021）：
```
v_imp = v_esc(M_total) = √(2G × M_total / R(M_total))
```
- late-stage GI は低速衝突が典型的（v_∞ ≈ 0）
- Canup (2004): 月形成 GI も v_imp ≈ v_esc
- Nakajima et al. (2021): N 体計算で v_imp ≤ 1.5 v_esc，中央値は v_esc 近辺
- R(M) は Sakuraba 2021 Eq.(3)（Seager et al. の bulk composition フィット）で算出

**温度補正された融解エネルギー E_m**（de Vries 2016 Eq.3, Abramov et al. 2012 より）：
```
E_m = E_m0 × [1 − C_p × (T_s + dT/dz × d_m) / (C_p × (T_l0 + dT_l/dP × P_dm) + L_m)]
```
- d_m : melt pool 深さ（数値反復で解く，後述）
- P_dm : 深さ d_m での圧力 = ρ_melt × g × d_m
- E_m は target 物質 1 kg を衝突後の減圧過程で溶融させるために必要な有効比エネルギー（J/kg）である。本モデルでは dunite 近似した target mantle に対する値を用いる。

**物性値**（de Vries 2016 Table 2; Pierazzo et al. 1997; Clauser 2011; Navrotsky 1995）：

| 変数 | 値 | 出典 |
|---|---|---|
| E_m0 | 9.0 × 10⁶ J/kg | dunite（Pierazzo et al. 1997） |
| C_p | 1300 J/kg/K | dunite melt（Clauser 2011） |
| L_m | 718 × 10³ J/kg | dunite（Navrotski 1995; de Vries 2016 Table 2） |
| T_s | 1750 K | 衝突前表面温度（de Vries 2016 Table 1） |
| T_l0 | 1950 K | 表面での液相線温度（de Vries 2016） |
| dT_l/dP | 28.3 K/GPa（P < 60 GPa），14.0 K/GPa（P ≥ 60 GPa） | Liebske & Frost 2012; de Vries 2016 Table 2 |
| dT/dz | 0.1 K/km | 下部マントル温度勾配（Monnereau & Yuen 2002; de Vries 2016） |
| ρ_mantle | ≈ 3300 kg/m³ | — |
| ρ_core | ≈ 8000 kg/m³ | — |

- Section 5 の質量収支では，de Vries の `V_melt` を **target 由来の silicate melt volume** と解釈する．impactor silicate はこの melt pool に全量取り込まれると近似し，impactor metal は別途金属 reservoir として加える

---

### 3-3. Melt Pool 深さ d_pool の計算

**melt pool を球形ジオメトリと仮定**（de Vries 2016; Pierazzo et al. 1997）：
```
V_melt = (π/3) × d_pool² × (3R_planet − d_pool)   ← 球冠体積公式（高さ d_pool，球半径 R_planet）
```
- 検算：d_pool = 2R → V = 4πR³/3（全球）✓，d_pool = R → V = 2πR³/3（半球）✓
**反復アルゴリズム（d_m と E_m の循環依存を解く）：**

E_m は d_m に依存し（Eq.3），V_melt は E_m に依存し（Eq.2），d_m は V_melt に依存する（球冠公式）→ 1変数 d_m の不動点反復で解く：

```
【初期値】
d_m^(0) = 0
   → E_m^(0) = E_m0 × [1 − C_p × T_s / (C_p × T_l0 + L_m)]   （深さ補正なしの値）

【反復】n = 0, 1, 2, ...
   (1) E_m^(n) = E_m0 × [1 − C_p×(T_s + dT/dz×d_m^(n)) / (C_p×(T_l0 + dT_l/dP×ρ_mantle×g×d_m^(n)) + L_m)]
   (2) V_melt^(n) = (π/6)×k×E_m^(−3μ/2)×(ρ_p/ρ_t)×D_p³×v_imp^(3μ)×sin^(2p_sin)(θ)
   (3) d_m^(n+1) = solve_spherical_cap(V_melt^(n), R_planet)
       （(π/3)×d²×(3R−d) = V を d について Newton 法で解く）

【収束判定】
   |d_m^(n+1) − d_m^(n)| < 1 m  → 収束
   最大反復回数：100 回（非収束時はエラー出力）
```

物理的フィードバック：`d_m ↑` に対し，液相線の深さ勾配 `dT_l/dz = (dT_l/dP)ρg` の方が地温勾配 `dT/dz` より大きい．したがって深くなるほど `T_liq(d_m)` の増加が `T_geotherm(d_m)` を上回り，融解に必要な比エネルギー `E_m` は増大する．Eq.2 より `V_melt ∝ E_m^(−3μ/2)` なので，`E_m ↑ → V_melt ↓ → d_m ↓` となり，これは負のフィードバックである．Earth-like な `ρg` を入れると `dT_l/dz ≈ 0.91 K/km`（`dT_l/dP = 28.3 K/GPa` 側）であり，`dT/dz = 0.1 K/km` より十分大きい．したがって固定点反復は通常 5–10 回で収束する．

- 上式と Eq.2（V_melt）を連立して d_pool を数値的に反復解（d_m を更新しながら E_m も更新）

**melt pool vs global MO の圧力の違い：Nakajima et al. (2021) の重要な知見：**

> *"The equilibrium pressure at the base of a melt pool can be higher (up to ≈ 80%) than those of radially-uniform global magma ocean models."*

- 等方拡散（isostatic adjustment）の timescale = 10²〜10⁵ 年（Reese & Solomatov 2006）
- 金属-ケイ酸塩平衡の timescale = 数時間〜数ヶ月（Dahl & Stevenson 2010）
- → 元素分配が終わるまで MO は melt pool のまま ⇒ melt pool 底圧力を使うべき

**採用する圧力・温度と equilibration の仮定：**
```
d_pool = min(d_pool_raw, R_planet − R_core)    ← mantle 厚さでクランプ（後述）
P_eq   = ρ_mantle × g(M_total) × d_pool        （melt pool 底での hydrostatic 圧力）
T_eq   = T_liq(P_eq) = T_l0 + (dT_l/dP) × P_eq （melt pool 底 = 液相線，Liebske & Frost 2012）
```
- global MO の平均深さではなく，melt pool 底を使う（Nakajima et al. 2021）
- 金属は底に向かって降下しながら逐次 re-equilibration を繰り返し，最終的に底での条件が支配的になる（Nakajima et al. 2021）

**P_eq の近似の正当化：**
- de Vries (2016) は melt pool 底での P を用いると述べるが，depth → pressure の変換式を明示していない
- 本モデルでは ρ_mantle = 3300 kg/m³ 一定とした単純 hydrostatic 近似（P = ρg d）を採用
- これは金属-ケイ酸塩平衡研究の標準的な慣行である（Rubie et al. 2015, *Icarus*; Nakajima et al. 2021）
- ρ の深さ依存性（下部マントル ~4400 kg/m³）を無視した誤差は P_eq の過小評価につながるが，感度解析で確認する

**T_eq = T_liq の正当化：**
- de Vries (2016) p.5：*"It has been assumed, in almost all studies of siderophile element partitioning, that the pressure of metal-silicate equilibration corresponds to the base of the magma ocean, with temperature being defined by the equivalent peridotite liquidus or solidus."*
- dT_l/dP の値は de Vries (2016) Table 2（Liebske & Frost 2012 より）で確認：28.3 K/GPa (P < 60 GPa)，14.0 K/GPa (P ≥ 60 GPa)，T_l0 = 1950 K
- **equilibration efficiency ξ_eq = 1（完全 re-equilibration）を仮定**
  - ξ_eq < 1 は物理的に正当化できるパラメータ化が難しいため採用しない（佐々木さんとの議論に基づく方針決定）

**melt pool がコアに届かない場合の metal 挙動と N/C/H の輸送機構：**

1. **melt pool 内での metal-silicate 平衡（軽元素の金属相への溶解）**
  
   impactor 由来の金属液滴（~1 cm，Dahl & Stevenson 2010）が melt pool 内を降下しながら silicate melt と平衡を取る．このとき N, C, H が金属相に化学的に溶解する：
   - N：鉄格子間原子として溶解（鉄ニトライド）
   - C：溶解 C または Fe₃C（鉄炭化物）として溶解
   - H：溶解 H として溶解
   
   D_N, D_C, D_H はまさにこの「melt pool 内での溶解」を定量する量である

2. **metal pond 形成とコアへの輸送**
   液滴が melt pool 底に集積 → metal pond 形成 → ダイアピル不安定によりコアへ沈降（timescale ~10³–10⁶ 年）→ N, C, H を金属相に乗せたままコアへ運搬

3. **固体マントル中での再平衡はない**
   > *"subsequently descend as large diapirs to the growing core without further equilibration with the surrounding silicate"*（de Vries 2016; Nakajima et al. 2021）
   
   固体との拡散交換は極めて遅く，降下中の再平衡は無視できる

- self-oxidation（Fe²⁺ → Fe³⁺ + Fe⁰）も同様に melt pool 内の高圧反応として成立：Fe³⁺ → ケイ酸塩メルトに残留，Fe⁰ → metal pond へ → ダイアピルでコアへ
- GI 間隔（~10⁷ 年）より timescale が十分短いため，mass balance 上は即時コア合体として扱う

**d_pool のクランプ（実装上の上限）：**
- d_pool > R_planet − R_core の場合は d_pool = R_planet − R_core にクランプ
- → 全 mantle が溶融した状態に相当し，P_eq = P_CMB（コアマントル境界圧力）となる

### 3-5. Fe³⁺/ΣFe(P, T, fO₂) の計算（Zhang 2024 SI Eq.S18 = Hirschmann 2022 Eq.21）

**採用モデルの根拠：**
- Armstrong et al. (2019): 10–23 GPa で disproportionation を初めて実験的に確認，ただし過大推定（Tait EOS 使用）
- Kuwahara et al. (2023): peridotite 組成で 15–28 GPa を実験，過大推定
- Zhang et al. (2024): 38–71 GPa の LH-DAC 実験（最高圧）で "fit 1" モデルを更新
  - Armstrong と Kuwahara は過大推定，Zhang 2024 が最も信頼できる
  - BSE 収支（Fe³⁺/ΣFe = 0.075–0.115）に整合するのも Zhang 2024 のみ

**基本方程式**（Zhang 2024 SI Eq.S18 = Hirschmann 2022 Eq.21）：
```
log(X_FeO₁.₅ / X_FeO) = a·log(fO₂) + b + c/T
  − ΔCp/(R·ln10) × [1 − T₀/T − ln(T/T₀)]
  − ∫_{P₀}^{P} ΔV dP / (RT·ln10)
  + (1/T) × [Y₁·X_SiO₂ + Y₂·X_TiO₂ + Y₃·X_MgO + Y₄·X_CaO
             + Y₅·X_NaO₀.₅ + Y₆·X_KO₀.₅ + Y₇·X_PO₂.₅
             + Y₈·X_SiO₂·X_AlO₁.₅ + Y₉·X_SiO₂·X_MgO]

r_ferric ≡ X_FeO₁.₅ / X_FeO
Fe³⁺/ΣFe = r_ferric / (1 + r_ferric)   （化学量論：FeO₁.₅ は 1 Fe³⁺/f.u.）
```

- `r_ferric` は Zhang/Hirschmann 式が直接与える ferric/ferrous 比であり，半径 `R(M)` とは無関係である

**self-oxidation に伴う Fe⁰ 生成の mass balance**

高圧シリケートメルト中の disproportionation は，単陽イオン酸化物基準では
```
3FeO(silicate) → Fe(metal) + 2FeO₁.₅(silicate)
```
で表される．したがって silicate 中の ferric iron が 2 mol 増えるごとに，1 mol の `Fe⁰` が metal pond / core へ移る．

反応進行度 `ξ_ox [mol]` を bookkeeping variable として導入すると，
```
n_FeO,new      = n_FeO,old − 3ξ_ox
n_FeO₁.₅,new   = n_FeO₁.₅,old + 2ξ_ox
N_Fe,sil,new   = N_Fe,sil,old − ξ_ox
```
ここで `N_Fe,sil = n_FeO + n_FeO₁.₅` は silicate melt 中に残っている total Fe のモル数である．`ξ_ox` は`Fe³⁺/ΣFe` の変化と上の化学量論から一意に決まる反応進行度である．

`f_old ≡ (Fe³⁺/ΣFe)_old`, `f_new ≡ (Fe³⁺/ΣFe)_new` とおくと，
```
ξ_ox        = N_Fe,sil,old × (f_new − f_old) / (2 + f_new)
M_Fe0,oxid  = ξ_ox × M_Fe
```
であり，生成した `Fe⁰` の質量 `M_Fe0,oxid` を追加コア質量として算入する．低圧では `f_new − f_old` が小さいため，`M_Fe0,oxid` も自然に 0 に近づく．

**Hirschmann 2022 Table 2 パラメータ：**

| パラメータ | 値 | 意味 |
|---|---|---|
| a | 0.1917 | fO₂ 指数 |
| b | 1.961 | 標準 ΔG 寄与 |
| c | 4158.1 K | 温度依存 ΔG |
| ΔCp | 33.25 J/(mol·K) | — |
| T₀ | 1673.15 K | 参照温度 |
| Y₁ (SiO₂) | 520.46 K | — |
| Y₂ (TiO₂) | 185.37 K | — |
| Y₃ (MgO) | 494.39 K | — |
| Y₄ (CaO) | 1838.34 K | — |
| Y₅ (NaO₀.₅) | 2888.48 K | — |
| Y₆ (KO₀.₅) | 3473.68 K | — |
| Y₇ (PO₂.₅) | 4473.6 K | — |
| Y₈ (SiO₂·AlO₁.₅) | 1245.09 K | — |
| Y₉ (SiO₂·MgO) | 1156.86 K | — |

**溶融体組成 X_i**（BSE 組成を単陽イオン基準モル分率に換算，McDonough & Sun 1995）：

| 成分 | X_i | 成分 | X_i |
|---|---|---|---|
| SiO₂ | 0.380 | AlO₁.₅ | 0.045 |
| MgO | 0.479 | CaO | 0.031 |
| FeO | 0.056 | NaO₀.₅ | 0.007 |
| TiO₂ | 0.001 | KO₀.₅, PO₂.₅ | ≈ 0 |

- FeO は self-oxidation 進行とともに逐次更新する．
- その他の major oxides（SiO₂, MgO, CaO, NaO₀.₅, TiO₂ など）は BSE ペリドタイトメルト代表値で固定する．self-oxidation の主反応は Fe²⁺ → Fe³⁺ + Fe⁰ であり，これらは Zhang/Hirschmann 式では組成補正項にのみ入る二次効果だからである
- 想定される FeO 変化量（Δ~0.007）に伴う他成分モル分率の再正規化効果は `~10^-3 log units` 程度と小さく，Fe³⁺/ΣFe の圧力依存や ΔIW_eq の不確実性より十分小さいため，現段階では固定近似でよい

**ΔV 圧力積分：Birch-Murnaghan EOS（Zhang 2024 SI Eq.S21）**

```
ΔV = V_FeO₁.₅(P,T) − V_FeO(P,T) を P₀ = 0.0001 GPa から P_eq まで数値積分
V_FeO₁.₅ ≡ V_Fe₂O₃/2  （化学量論より）
```

- Fe₂O₃ と FeO の圧縮率の違いにより，高圧ほど ΔV が負になり Fe³⁺ が安定化される
- この項こそが「圧力による self-oxidation」の物理的本体であり，省略すると Fe³⁺/ΣFe の圧力依存性が消えて self-oxidation が再現できない

各酸化物の体積 V(P,T) は 4 次 BM-EOS + 熱圧力補正で計算：
```
P(V,T)     = P(V,T_ref) + B_th(V) × (T − T_ref)   （Zhang 2024 SI Eq.S20）
P(V,T_ref) ← 4 次 BM-EOS                           （Zhang 2024 SI Eq.S21）
T_ref      = 3000 K

B_th(V) = [a − b×(V/V₀) + c×(V/V₀)²] / 1000   [GPa/K]   （Zhang 2024 SI Eq.S22; Deng et al. 2020）
```

**Fit 1 パラメータ**（Zhang 2024 SI Table S4 "fit 1"）：

| 相 | V₀ [Å³] | K₀ [GPa] | K' | K'' [GPa⁻¹] |
|---|---|---|---|---|
| Fe₂O₃（per f.u.） | 1204.69 | 21.35(14) | 3.626(29) | 0.009（fixed） |
| FeO（per f.u.） | 1180.10 | 24.22(16) | 3.382(32) | 0.012（fixed） |

V₀ と K'' は Deng et al. (2020) 12.5 mol% FeO モデルから固定，K₀ と K' のみフィット．

**熱圧力係数 a, b, c**（Deng et al. 2020 12.5% モデルから固定；B_th [GPa/K] × 1000 の単位）：

| 相 | a | b | c |
|---|---|---|---|
| Fe₂O₃ | 34.53 | 68.64 | 35.27 |
| FeO | 35.70 | 71.10 | 36.60 |

実装：V(P,T) は P→V の逆関数を Newton 法で数値的に求め，ΔV(P) を台形則 or Gauss 求積で P₀ → P_eq まで積分．

---

### 3-6. ΔIW の計算：2 種類の使い分け

> **重要な区別**：(A) 金属-ケイ酸塩平衡中の ΔIW（分配係数用）と，(B) 金属沈降後の表面 ΔIW（大気化学用）は異なる量であり，別々に扱う．

#### (A) ΔIW_eq,deep：melt pool 底での金属-ケイ酸塩平衡条件（D_N, D_H, D_C の入力）

反応（Zhang 2024 SI Eq.S23）：
```
FeO(silicate liquid) → Fe(molten alloy) + 1/2 O₂
```

活量（Zhang 2024 SI Eqs.S24–S25; Hirschmann 2022）：
```
a_FeO^silicate = γ_FeO × X_FeO
  γ_FeO = 1.55  （Zhang 2024 SI Eq.S25；Hirschmann 2022 の高圧実験コンパイルに基づく代表値）
a_Fe^alloy = 0.75  （Zhang 2024 SI Eq.S25；PtFe 合金実験値の代表値）
```

Zhang 2024 の数値結果（Komabayashi 2014 EOS を使った計算結果）では，high-T / low-T adiabat の両方で `ΔIW_eq` は `~−2 ± 0.15` に収まる．したがって **初期値** として `ΔIW_eq,deep = −2.0` を採用するのは妥当である．ただし本モデルでは，これを全 step の固定値にはせず，**各衝突ごとに self-oxidation 後の melt pool 底状態から更新する**．

各衝突 step では，melt pool 底 (`P_eq, T_eq`) において次の fixed-point を解く：

1. まず前 step から引き継いだ `ΔIW_eq,deep^(n)` を仮定する（初回のみ `−2.0`）
2. Zhang 2024 / Hirschmann 2022 の `Fe³⁺/ΣFe(P,T,fO₂)` 式から，その `ΔIW_eq,deep^(n)` に対応する平衡 ferric fraction を計算する
3. `3FeO → Fe + 2FeO₁.₅` の化学量論で self-oxidation を適用し，`M_Fe0,oxid`，`X_FeO`，`Fe³⁺/ΣFe` を更新する
4. 更新後の ferrous FeO 活量 `a_FeO^silicate = γ_FeO X_FeO` から，`ΔIW_eq,deep^(n+1) = 2log10(a_FeO^silicate / a_Fe^alloy)` を計算する
5. `|ΔIW_eq,deep^(n+1) − ΔIW_eq,deep^(n)| / max(1, |ΔIW_eq,deep^(n)|) < 10^-6` で収束とする

すなわち，本モデルで `D_N`, `D_H`, `D_C`, `Δ¹⁵N^(metal-silicate)` に入るのは，**固定 `−2` ではなく，その衝突で self-oxidation と整合した melt pool 底の `ΔIW_eq,deep`** である．

実装上の順序は
```
inherited ΔIW_eq,deep
→ Fe³⁺/ΣFe equilibrium at (P_eq, T_eq)
→ self-oxidation update
→ activity-based update for ΔIW_eq,deep
→ D_N, D_H, D_C, Δ¹⁵N^(metal-silicate)
```
である．この `ΔIW_eq,deep` は次 step の初期 deep redox state として planet mantle に保持する．

#### (B) ΔIW_surface：金属沈降後の表面 fO₂（大気化学・outgas用）

Zhang et al. (2024, SI) Eq.S18 を P ≈ 0，T = T_surface で評価し，Fe³⁺/ΣFe から逆算：
```
log(fO₂)      = [log(X_FeO₁.₅/X_FeO) − b − c/T − ΔCp 項 − 組成項] / a
ΔIW_surface   = log(fO₂) − log(fO₂_IW(T_surface, P ≈ 0))
```

- ここで求める `ΔIW_surface` は，表面で metal-silicate が新たに平衡して決まる redox condition ではなく，**深部 melt pool で確立した ferric ratio をもつ melt が表面に到達したときに，その組成が見かけ上どの fO₂ に対応するか**を表す診断量である
- したがって `ΔIW_surface` は `D_N`, `D_H`, `D_C` に入れる平衡入力ではなく，大気主成分の選択や outgassing redox の記述に用いる二次的な量である
- この逆算では，金属沈降後の melt が上昇中に直ちに bulk redox を再平衡すると仮定せず，`P_eq, T_eq` で成立した `Fe³⁺/ΣFe` を保ったまま表面へ運ばれると近似する（**kinetic quench 仮定**）
- これはモデル近似であり，対流上昇中の部分的な再平衡は現段階では明示的に扱わない．必要なら感度テストで評価する

**T_surface の計算（断熱プロファイルからの外挿）：**

melt pool 内部は対流により断熱的であり，底（T_eq, P_eq）から表面（P = 0）まで断熱線に沿って温度が下がる：
```
T_surface = T_eq − (dT_a/dP) × P_eq
```

- **dT_a/dP = 15 K/GPa**（シリケートメルトの断熱温度勾配代表値；Wolf & Bower 2018, EPSL; Stixrude et al. 2009, EPSL）
  - 上部マントルメルト：~20–30 K/GPa，下部マントルメルト：~10–15 K/GPa → 一定値 15 K/GPa を全圧力範囲の代表値として採用
- 下限クランプ：T_surface = max(T_eq − (dT_a/dP) × P_eq, T_l0)
  - melt pool 表面が液相線（T_l0 = 1950 K）を下回ることを防ぐ

#### (C) IW buffer の絶対 fO₂：Hirschmann (2021) Table 1

定義より：
```
log(fO₂) = ΔIW + log(fO₂_IW(T, P))
```

log(fO₂_IW) の経験式（Hirschmann 2021 Table 1；T [K]，P [GPa]）：
```
fcc/bcc 鉄が安定な領域：  log fO₂ = a + b·T + c·T·ln(T) + d/T
hcp 鉄が安定な領域：      log fO₂ = e + f·T + g·T·ln(T) + h/T
```

各係数の P 依存性（m = m₀ + m₁P + m₂P² + m₃P³ + m₄P^(1/2)）：

| 係数 | m₀ | m₁ | m₂ | m₃ | m₄ |
|---|---|---|---|---|---|
| a | 6.844864 | 1.175691E-1 | 1.143873E-3 | 0 | 0 |
| b | 5.791364E-4 | -2.891434E-4 | -2.737171E-7 | 0 | 0 |
| c | -7.971469E-5 | 3.198005E-5 | 0 | 1.059554E-10 | 2.014461E-7 |
| d | -2.769002E+4 | 5.285977E+2 | -2.919275E+0 | 0 | 0 |
| e | 8.463095 | -3.000307E-3 | 7.213445E-5 | 0 | 0 |
| f | 1.148738E-3 | -9.352312E-5 | 5.161592E-7 | 0 | 0 |
| g | -7.448624E-4 | -6.329325E-6 | 0 | -1.407339E-10 | 1.830014E-4 |
| h | -2.782082E+4 | 5.285977E+2 | -8.473231E-1 | 0 | 0 |

**実装上の場合分け**（Fe 結晶相の物理的意味は意識せず，圧力で式を切り替えるだけ）：
```
P_boundary [GPa] = −18.64 + 0.04359·T − 5.069×10⁻⁶·T²
P > P_boundary  → hcp 式を使用
P ≤ P_boundary  → fcc/bcc 式を使用
```
- 適用範囲：1000–3000 K，0.0001–100 GPa
- 精度：r.m.s. = 0.0065 log units，最大誤差 0.028 log units

---

## 4. 元素および同位体の分配

### 4-1. N（Shi et al. 2022）

**分配係数 D_N（Shi 2022 Eq.1）：**
```
log(D_N^(metal/silicate)) = −0.38 + 3370/T + 75×P/T
  + 0.56×NBO/T + 0.47×ΔIW
  + 1.97×log(1 − x_S^metal) + 5.73×log(1 − x_C^metal)
  + 3.88×log(1 − x_Ni^metal) + 5.52×log(1 − x_Si^metal)
```

- T [K]，P [GPa]
- **ΔIW = self-oxidation と activity 比から更新された `ΔIW_eq,deep`**（Section 3-6A）
- **NBO/T ≈ 2.5**：BSE 組成（McDonough & Sun 1995）から以下の標準式で計算：
  ```
  NBO = 2n_Mg + 2n_Fe²⁺ + 2n_Ca + n_Na + n_K − 2n_Al   （modifier 寄与 − Al 消費分）
  T   = n_Si + n_Al + n_Ti                               （tetrahedral cation 数）
  NBO/T = NBO / T

  BSE 値（per 100 g）：n_Si=0.747, n_Mg=0.938, n_Fe=0.106, n_Al=0.087,
                       n_Ca=0.063, n_Na=0.012, n_Ti=0.005
  → NBO = 2.052, T = 0.839, NBO/T ≈ 2.45 ≈ 2.5
  ```
  Shi 2022 の実験範囲 0.02–3.12 のほぼ中央．各ステップで X_i を更新する場合は再計算可能だが，FeO の変化量（Δ~0.007）に対する NBO/T の感度は低いため固定値を採用
- `x_i^metal` は，その衝突で silicate melt と平衡する **impactor 金属相の alloy composition** `x_i^imp_met` を表す．これは地球 core 全体の平均組成とは異なる

金属メルト組成の扱い：
- `x_S^imp_met`, `x_Ni^imp_met`, `x_Si^imp_met` は分化済み EC impactor の representative alloy composition として固定する
- `x_C^imp_met` は 0.01→0.10 `M_Earth` の precompute における**初期値のみ**感度パラメータとし，以後はその precompute の出力コア組成から更新する
- GI1–GI7, GI9, GI8, late veneer では，同じ 0.10 `M_Earth` 時点の precompute 出力 alloy を proxy として用いる

初期値（分化済み EC impactor representative alloy composition）：

| 成分 | 初期値 |
|---|---|
| x_Si^metal | ≈ 0.06（Si 核移送の EC モデル代表値） |
| x_Ni^metal | ≈ 0.05 |
| x_S^metal | ≈ 0.02 |
| x_C^metal | **precompute 初期値のみ感度パラメータ**（デフォルト 0.01；範囲 0.005–0.05 で感度テスト） |

x_C^metal の設定理由：EC 金属相の炭素は主に cohenite (Fe₃C) や graphite として独立固相で存在し，Fe-Ni 合金への溶存量の直接測定値が文献にない．そのため**0.01→0.10 precompute の開始時点**に与える `x_C^imp_met` の初期値のみ感度パラメータとして扱う．デフォルト値 0.01 は EC バルク C（~0.3–0.5 wt%）と Fe-C 二元系溶解度（Dasgupta & Hirschmann 2010）から保守的に設定した推定値である．この初期値を用いた precompute の出力コア組成から `x_C^imp_met(0.10 M_Earth)` を求め，GI1–GI7 および GI9 の `D_N`, `D_H` 評価に用いる（Section 5）．

GI8（CC event）についても，`D_N`, `D_H`, `D_C` に入る金属相 alloy 組成は EC-like proxy を用いる．具体的には `x_S`, `x_Ni`, `x_Si` は同じ representative values を保ち，`x_C` は 0.10 `M_Earth` precompute 出力コア組成から得た値を用いる．late veneer についても，同じ precompute-derived alloy proxy を用いる．ここで bulk volatile inventory の CI-like 仮定と alloy proxy は役割が異なる量であり，後者の近似誤差は感度テストで評価する．

**N 同位体分別係数（Shi 2022 Eq.2）：**
```
Δ¹⁵N^(metal-silicate) = (−1.6×10⁷ − 1.1×ΔIW×10⁷) / T²
```
- Δ¹⁵N^(metal-silicate) = δ¹⁵N^metal − δ¹⁵N^silicate [‰]
- 圧力・NBO/T・軽元素含有量の影響は統計的に無視可能（Shi 2022 本文）

マントル δ¹⁵N への反映：各衝突ステップで，コアへ移る金属メルトの δ¹⁵N = δ¹⁵N^silicate + Δ¹⁵N^(metal-silicate)．残留マントルの δ¹⁵N は質量収支で更新する．

---

### 4-2. H（Tsutsumi et al. 2025）

**分配係数 D_H（Tsutsumi 2025 Eq.3）：**
```
log₁₀(D_H^(metal/silicate)) + ΔIW/4 = a + b/T + c×P/T + d×log₁₀(1 − x_C^metal)
```
整理すると：
```
log₁₀(D_H) = a + b/T + c×P/T + d×log₁₀(1 − x_C^metal) − ΔIW/4
```

フィット定数（Tsutsumi 2025 フィット結果）：

| 定数 | 値 |
|---|---|
| a | 2.31(20) |
| b | −1542(500) |
| c | −38.9(78) |
| d | 6.69(30) |

- T [K]，P [GPa]，x_C^metal：金属中 C のモル分率，ΔIW：Section 3-6A の値
- **炭素による非理想相互作用（d 項）が重要**：x_C^metal が大きいほど D_H が増加
- この式は Tagawa et al. の C-poor 系データと本研究の C-rich 系データを同時フィット

---

### 4-3. C（Fischer et al. 2020）

**分配係数 D_C（Fischer 2020 Eq.1，nanoSIMS 版）：**
```
log₁₀(D_C) = 1.49 + 3000/T − 235×P/T
  + 9.6×log₁₀(1 − X_S) − 19.5×log₁₀(1 − X_O)
  − 0.118×NBO/T − 0.238×ΔIW
```

- T [K]，P [GPa]
- X_S, X_O：金属メルト中の S および O のモル分率
- **ΔIW = self-oxidation と activity 比から更新された `ΔIW_eq,deep`**（Section 3-6A）
- **NBO/T ≈ 2.5**（Section 4-1 と共通）

**X_S の扱い：**
- `X_S` も equilibrating impactor metal の alloy composition を表す量であり，ここでは `X_S ≈ 0.02` の representative value を用いる

**X_O の計算（Fischer et al. 2015 Table 2 の O exchange coefficient を使用）：**

X_O は独立した溶解度式ではなく，Fischer et al. (2015, *GCA*) が Si/O 分配に対して与えた O の exchange coefficient `K_D^O` から計算する：

```
K_D^O = (X_O^metal × X_Fe^metal) / X_FeO^silicate
log10(K_D^O) = a_O + b_O/T + c_O×P/T
```

Fischer et al. (2015) Table 2 では，O に対して
```
a_O = 0.6,
 b_O = −3800 K,
 c_O ≈ 0
```
であり，O partitioning は `P > 5 GPa` のデータ範囲で圧力依存が弱い．したがって本モデルでは

```
log10(K_D^O) = 0.6 − 3800/T
X_O^metal    = K_D^O × X_FeO^silicate / X_Fe^metal
```

を用いる．

- **X_FeO^silicate**：Section 3-5 で self-oxidation を通じて逐次更新される量（初期値 0.056）
- **X_Fe^metal** = 1 − x_S − x_C − x_Ni − x_Si（equilibrating impactor metal の Fe モル分率）
- この形は Fischer 2015 が直接フィットした `K_D^O` を用いるものであり，`X_O` を別の経験式から外挿するよりも本モデルの圧力温度範囲に対して保守的である

重要な注意：
- Fischer 2020 Eq.(1) は C-undersaturated 条件（graphite なし）のデータに基づく → 本モデルの想定（EC-like impactor で C 供給があるが過飽和でない）に整合
- `D_C` に入る `X_O` は上の `K_D^O` から毎ステップ自己整合的に更新する

---

### 4-4. 共通入力変数のまとめ

| 変数 | 値 / 更新方法 | 出典 |
|---|---|---|
| T | 各ステップの melt pool 温度（T_liq(P_eq)） | de Vries / Liebske & Frost 2012 |
| P | P_eq（Section 3-3） | — |
| ΔIW | self-oxidation と activity 比から更新された ΔIW_eq,deep | Section 3-6A |
| NBO/T | ≈ 2.5（固定；ペリドタイトメルト代表値） | Shi 2022 |
| x_S^metal | equilibrating impactor metal の representative value（固定 ≈ 0.02） | — |
| x_C^metal | precompute 初期値のみ感度パラメータ；以後は 0.01→0.10 出力コア組成から決定 | — |
| x_Ni^metal | equilibrating impactor metal の representative value（固定 ≈ 0.05） | — |
| x_Si^metal | equilibrating impactor metal の representative value（固定 ≈ 0.06） | — |
| X_S (Fischer) | equilibrating impactor metal の representative value（固定 ≈ 0.02） | — |
| X_O (Fischer) | 平衡から毎ステップ計算：K(P,T)×X_FeO/X_Fe | Fischer et al. 2015 |

---

## 5. 質量収支方程式

各衝突ステップ（impactor 質量 M_imp，鉄コア比率 f_Fe）での更新手順．

**f_Fe の定義と採用値**

- `f_Fe` は，そのステップで金属相として分離・沈降する **impactor 由来 metallic Fe mass fraction** を表す
- 分化済み impactor では core mass fraction，未分化 impactor では bulk 中にもともと含まれる金属 Fe fraction として解釈する
- デフォルトでは Sakuraba.C の `X_m = 0.325` に合わせ，`f_Fe = 0.325` を全フェーズの representative value とする
- GI8（CC embryo）についても direct constraints が弱いため，まずは同じ `f_Fe = 0.325` を proxy とし，必要なら感度テストで振る

**単位系の約束**

- 貯留層濃度 `C_mantle,i`, `C_core,i`, `C_imp,i`, `C_imp_met,i`, `C_sil,i`, `C_melt,i` はすべて **ppm（質量濃度）** とする
- 保存式の中では `w_i ≡ C_i × 10^-6` として質量分率に変換して扱う
- target 由来 melt mass は `M_from_tgt = ρ_melt × V_melt` で与え，`ρ_melt ≈ ρ_mantle = 3300 kg/m³` を用いる
- `N_A = 6.02214076 × 10^23 mol^-1`（Avogadro 定数）
- `M_t = M_core_old + M_mantle_old` を current target mass，`R_planet = R(M_t)` をその半径とする
- `P_eq` は GPa，`T` は K，atmospheric partial pressure `P_i` は Pa とし，Henry 型溶解度式では `(P_i / 1 MPa)` を用いて無次元化する
- atmosphere reservoir は元素質量 `M_atm,C`, `M_atm,H`, `M_atm,N` で追跡するが，partial pressure 計算では各 gas species の分子量へ換算する
- 各 gas species について，元素質量 `M_atm,i` から分子数 `N_i` へは `N_i = M_atm,i / m_elem,i` で変換する
  - `m_elem,C = 12/N_A`（CO₂, CO, CH₄ のいずれでも 1 molecule あたり C は 12 amu）
  - `m_elem,H = 2/N_A`（H₂O, H₂ のいずれでも 1 molecule あたり H は 2 amu）
  - `m_elem,N = 28/N_A`（N₂ 1 molecule あたり N は 28 amu）
- atmosphere の mean molecular mass は `m̄ = M_atm,total / (N_C + N_H + N_N)` とする
  - `M_atm,total = ν_C M_atm,C + ν_H M_atm,H + ν_N M_atm,N`
  - `ν_C = 44/12, 28/12, 16/12` for CO₂, CO, CH₄
  - `ν_H = 18/2, 1` for H₂O, H₂
  - `ν_N = 1` for N₂

**貯留層の定義：**

| 貯留層 | 変数 |
|---|---|
| mantle silicate | M_mantle，C_mantle,i [ppm]，δ¹⁵N_mantle |
| core | M_core，C_core,i [ppm]，δ¹⁵N_core |
| atmosphere | M_atm,i [kg, 元素質量]，δ¹⁵N_atm |

### Step 1：Mixing（melt pool 形成）

```
M_from_tgt = min(ρ_melt × V_melt, M_mantle_old)          （溶融した target mantle）
M_imp_met  = f_Fe × M_imp                                （impactor metallic Fe）
M_imp_sil  = (1 − f_Fe) × M_imp                          （impactor silicate）
M_melt_sil = M_from_tgt + M_imp_sil                      （equilibrating silicate melt mass）
```

impactor の揮発性 inventory には 2 つのモードを区別する。

1. **phase-resolved impactor**（GI1–GI7, GI9）
```
w_melt,i⁰ = [w_mantle,i × M_from_tgt + w_imp,sil,i × M_imp_sil] / M_melt_sil
C_melt,i⁰ = 10^6 × w_melt,i⁰

N_melt     = w_melt,N⁰ × M_melt_sil
δ¹⁵N_melt = [w_mantle,N × M_from_tgt × δ¹⁵N_mantle
              + w_imp,sil,N × M_imp_sil × δ¹⁵N_imp] / N_melt
```

2. **bulk-inventory impactor**（暴走・寡占成長の未分化 EC，GI8，late veneer）
```
M_imp,vol,i = w_imp,bulk,i × M_imp
w_melt,i⁰   = [w_mantle,i × M_from_tgt + M_imp,vol,i] / M_melt_sil
C_melt,i⁰   = 10^6 × w_melt,i⁰

N_melt      = w_melt,N⁰ × M_melt_sil
δ¹⁵N_melt  = [w_mantle,N × M_from_tgt × δ¹⁵N_mantle
               + M_imp,vol,N × δ¹⁵N_imp] / N_melt
```

ここで `w_mantle,i = C_mantle,i × 10^-6`, `w_imp,sil,i = C_imp,sil,i × 10^-6`, `w_imp,bulk,i = C_imp,bulk,i × 10^-6` である．bulk-inventory impactor では，揮発性元素総量を `bulk inventory × M_imp` として保持したうえで，その全量を equilibrating silicate reservoir に投入する bookkeeping を採用する．これは未分化体の金属 Fe 質量が 0 であることを意味しない．金属 Fe 自体は `M_imp_met = f_Fe × M_imp` として別途存在し，melt pool 内で分離・沈降する。

- `V_melt` は target 側に新たに生成した silicate melt volume と解釈し，impactor silicate はこの melt pool に全量混合すると近似する
- したがって Step 1 で平衡に参加する silicate は `M_from_tgt + M_imp_sil`，metal は `M_imp_met` である

### Step 2：三貯留層 Equilibrium（Sakuraba 2021 Eq.7–9 を melt pool に適用）

melt pool 内で atm / silicate melt / metal が同時平衡を取る．

**保存式（平衡前 = 平衡後）：**
```
M_atm,i^old + w_melt,i⁰ × M_melt_sil + w_imp_met,i × M_imp_met
    = M_atm,i^new + w_sil,i × M_melt_sil + D_i × w_sil,i × M_imp_met
```

- 左辺：平衡前の総揮発性元素量（既存大気 + melt pool シリケート + impactor 金属相）
- 右辺：平衡後の三貯留層への分配（Sakuraba 2021 Eq.7 のアナログ）
- `D_i × w_sil,i`：金属-ケイ酸塩平衡後の金属の質量分率（Sakuraba 2021 Eq.8）
- ここで用いる `D_N`, `D_H`, `D_C` は，Section 3-5 / 3-6A の self-oxidation solve で収束した `ΔIW_eq,deep` を入力として計算する

**impactor 金属の揮発性元素濃度と同位体組成：**

| フェーズ | C_imp_met,i | δ¹⁵N_imp_met |
|---|---|---|
| 暴走・寡占成長（0.01→0.10） | 0（未分化体，pre-existing core なし） | 未使用 |
| GI1–GI7, GI9 = Theia | 0.01→0.10 計算の出力コア濃度 C_core,i | 0.01→0.10 計算の出力コア δ¹⁵N |
| GI8（CI-like CC embryo） | 0（bulk inventory は CI-like として silicate reservoir に一括投入し，pre-existing metal volatile reservoir は与えない） | 未使用 |
| late veneer | 0（未分化体として扱う） | 未使用 |

ここで `C_imp_met,i = 0` は「pre-existing な分化コアの揮発性濃度を持たない」ことを意味する．`M_imp_met = f_Fe × M_imp` 自体は，未分化 EC impactor や GI8 impactor に含まれる金属 Fe が melt pool 内で分離・沈降して形成する原始コア成分として別途存在しうる．`GI8` では，分配係数 `D_N`, `D_H`, `D_C` に入る alloy 組成 `x_S`, `x_Ni`, `x_Si`, `x_C` のみ EC proxy を用い，`w_imp_met,i = C_imp_met,i × 10^-6 = 0` とする．

**大気-ケイ酸塩平衡（Sakuraba 2021 Eq.9，Henry 則）：**
```
C_sil,i [ppm] = S_i × (P_i / 1 MPa)^(1/x_i)
w_sil,i       = C_sil,i × 10^-6
N_i           = M_atm,i^new / m_elem,i
m̄            = M_atm,total / (N_C + N_H + N_N)
P_i [Pa]      = g × m̄ × N_i / (4πR_planet²)
```

ここで `M_atm,total = ν_C M_atm,C + ν_H M_atm,H + ν_N M_atm,N` は atmosphere の総**分子**質量である．

保存式と Henry 則を連立すると，C と N については未知数 `C_sil,i`（または `P_i`）に対する **1 変数方程式**が得られる．ただし `m̄` が `M_atm,C`, `M_atm,H`, `M_atm,N` の全てに依存するため，Step 2 全体は以下の外側反復で self-consistent に解く：

1. 前ステップの `m̄` を初期値として与える
2. C, N は Henry 則＋保存式を 1 変数 root solve で解く
3. H は下記の H₂O/H₂ モデルで解く
4. 更新された `M_atm,C`, `M_atm,H`, `M_atm,N` から `m̄_new` を再計算する
5. `|m̄_new − m̄_old| / m̄_old < 10^-4` で収束，通常 3–5 回で十分

**gas species と溶解度 S_i**（ΔIW_surface からステップ毎に決定，Section 3-6B）：

| ΔIW_surface | C 種 | H 種 | N 種 |
|---|---|---|---|
| > −1 | CO₂ | H₂O | N₂ |
| −3 〜 −1 | CO | H₂ | N₂ |
| ≤ −3 | CH₄ | H₂ | N₂ |

Sakuraba et al. (2021) Table 1 の `S_i` は，上式のように `(P_i / 1 MPa)` を用いたとき `C_sil,i` を ppm で与える定数である．したがって本モデルでは `C_sil,i` を ppm のまま保持し，保存式の中だけ `w_sil,i = C_sil,i × 10^-6` に変換する．

- C は Sakuraba et al. (2021) Table 1 の nominal values `S_C = 1.6, 0.55, 0.22 ppm` をそれぞれ `CO₂, CO, CH₄` に対応させて用いる
- N は同じく Table 1 の nominal values `S_N = 1.0, 5.0, 50.0 ppm` を oxidized / intermediate / reduced regime に対応させて用いる
- N 種は現バージョンでは全 redox 条件で `N₂` に固定する．`NH₃`/HCN chemistry は範囲外とする
- H は Sakuraba.C と同様に，酸化条件で H₂O が主成分のときは Moore (1998) 型溶解度モデルを優先する．H 種が `H₂` のときは，Table 1 reduced-model row の published H₂ solubility `S_H = 5.0 (P_i/1 MPa)^(1/2) ppm` を Henry 型 `S_i` として用いる

**H₂O 条件での H 溶解度（Moore 1998；Sakuraba.C 実装準拠）**

H 種が `H₂O` のときは単純 Henry 則を使わず，silicate 中の H 質量 `M_H,sil` を未知数として次の方程式を root solve する：
```
0 = c_M log(P_H2O / 1 MPa)
    + (b_Al X_Al + b_Fe X_Fe + b_Na X_Na) (P_H2O + P_C + P_N) / (T × 1 MPa)
    − 2 log[(m_H2O × M_H,sil) / (r_CtoX × M_melt_sil × 2m_H)]
    + a_M/T + d_M
```

- 定数は Sakuraba.C と同じく `a_M = 2565`, `b_Al = −1.997`, `b_Fe = −0.9275`, `b_Na = 2.736`, `X_Al = 0.137`, `X_Fe = 0.124`, `X_Na = 0.0268`, `c_M = 1.171`, `d_M = −14.21`, `r_CtoX = 0.309`
- `m_H2O = 18/N_A`, `m_H = 1/N_A` を用い，圧力は MPa に換算して代入する
- `M_H,total = M_atm,H^old + w_melt,H⁰ × M_melt_sil + w_imp_met,H × M_imp_met`
- `P_H2O`, `P_C`, `P_N` は同じ外側反復で更新される atmospheric partial pressures
- 解が得られたら `M_H,met = D_H × M_H,sil × (M_imp_met / M_melt_sil)`，`M_H,atm = M_H,total − M_H,sil − M_H,met` とする
- H 種が `H₂` のときは C, N と同様に Henry 型 `S_i` を用いる

### Step 3：貯留層更新

**Core：**
```
M_core_new   = M_core_old + M_imp_met + M_Fe0,oxid
M_mantle_new = M_mantle_old + M_imp_sil − M_Fe0,oxid

w_core,i,new ← (w_core,i,old × M_core_old + D_i × w_sil,i × M_imp_met) / M_core_new
C_core,i,new ← 10^6 × w_core,i,new

M_core   ← M_core_new
M_mantle ← M_mantle_new
```

- `M_imp_met` は impactor にもともと含まれる金属 Fe の分離・沈降分，`M_Fe0,oxid` は self-oxidation により新たに生成した純鉄分である
- `M_Fe0,oxid` は N/C/H を運ばないため，`C_core,i` の分子には入らず，分母を通じて core 中揮発性濃度を希釈する
- `M_mantle_new` では，accreted silicate `M_imp_sil` を加え，silicate から core へ抜けた `M_Fe0,oxid` を差し引く

**0.01→0.10 precompute における impactor metal C モル分率の更新：**
更新後の embryo core C 濃度から，precompute の次ステップで用いる `x_C^imp_met` を再構築する：
```
w_C^core  = C_core,C × 10^-6                      （ppm → 質量分率）
w_Fe^core = 1 − w_C^core − w_S^core − w_Ni^core − w_Si^core

x_C^imp_met = (w_C^core / M_C)
              / [(w_Fe^core / M_Fe) + (w_C^core / M_C)
                 + (w_S^core / M_S) + (w_Ni^core / M_Ni) + (w_Si^core / M_Si)]
```

- 原子量：`M_Fe = 55.845`, `M_C = 12.011`, `M_S = 32.06`, `M_Ni = 58.693`, `M_Si = 28.085`
- `w_S^core, w_Ni^core, w_Si^core` を明示的に追跡する run ではその時点の値を用いる．簡略版では current representative alloy composition（`x_S^imp_met`, `x_Ni^imp_met`, `x_Si^imp_met`）に対応する値で近似してよい
- この再計算した `x_C^imp_met` は **0.01→0.10 precompute の次ステップ**の `D_N`（Shi 2022）と `D_H`（Tsutsumi 2025）に渡す
- 本計算の GI フェーズでは，地球 core の `C_core,C` を後続 impactor の `x_C^metal` に直接フィードバックしない．GI1–GI7, GI9 には `0.10 M_Earth` 時点の precompute 出力 `x_C^imp_met` を用いる

**Mantle**（N/H/C とも melt pool 分は C_sil,i を保持，well-mixed 平均で更新）：
```
w_mantle,i,new = (w_sil,i × M_melt_sil + w_mantle,i,old × (M_mantle_old − M_from_tgt)) / M_mantle_new
C_mantle,i_new = 10^6 × w_mantle,i,new

δ¹⁵N_mantle_new = (w_sil,N × M_melt_sil × δ¹⁵N_sil
                   + w_mantle,N,old × (M_mantle_old − M_from_tgt) × δ¹⁵N_mantle_old)
                  / (w_mantle,N,new × M_mantle_new)
```

- **Sakuraba 整合**：Henry 則で溶解している分（C_sil,i）はマントルに保持し，大気との分配は Step 2 の平衡結果（M_atm,N^new）に反映済み
- **「全量 outgas」は不採用**：GI 規模の melt pool（深さ数百〜数千 km）が固化した後，固体シリケート中の N 拡散係数は極めて低く（D_solid << D_melt），GI 間隔（~10⁷ 年）で N が抜けきらない → C_sil,N はマントルに保持される
- **well-mixed 仮定**：各 GI が異なる場所を溶かすため，melt pool は毎回マントルのランダムな部分をサンプリングするはず → 期待値として全体の well-mixed 平均と等価
- **M_melt_sil → M_mantle の極限**（大型 GI，全球溶融）：δ¹⁵N_mantle → δ¹⁵N_sil（Sakuraba の全球 MO モデルに一致）
- **M_melt_sil → 0 の極限**（小型衝突）：δ¹⁵N_mantle ≈ 変化なし

**Atmosphere**（Step 2 の平衡結果をそのまま使用，固化時の追加 outgas なし）：
```
M_atm,i = M_atm,i^new   （i = C, H, N；Step 2 の保存式の結果）
```

### Step 4：δ¹⁵N の同位体 mass balance

Shi et al. (2022) の3貯留層同位体収支（Eq.16, 19）を適用する．

**Δ¹⁵N^atm-silicate ≈ 0（Shi 2022 Eq.19）：**
高温での大気-シリケートメルト間の N 同位体分別は無視できる → **平衡後 δ¹⁵N^atm = δ¹⁵N^sil**

**δ¹⁵N_bulk の計算（全 N の質量加重平均）：**
```
Δ¹⁵N      = (−1.6×10⁷ − 1.1×ΔIW×10⁷) / T_eq²                          （Shi 2022 Eq.2）

N_total    = N_melt + M_atm,N^old + w_imp_met,N × M_imp_met

δ¹⁵N_bulk = (N_melt × δ¹⁵N_melt
              + M_atm,N^old × δ¹⁵N_atm^old
              + w_imp_met,N × M_imp_met × δ¹⁵N_imp_met) / N_total
```

**3貯留層同位体収支（Shi 2022 Eq.16，Δ¹⁵N^atm-sil = 0 を代入）：**
```
N_total × δ¹⁵N_bulk = N_total × δ¹⁵N_sil + N_met × Δ¹⁵N

→ δ¹⁵N_sil = δ¹⁵N_bulk − (N_met / N_total) × Δ¹⁵N
→ δ¹⁵N_met = δ¹⁵N_sil + Δ¹⁵N
→ δ¹⁵N_atm^new = δ¹⁵N_sil   （Δ¹⁵N^atm-sil = 0 より）
```

ここで
```
N_met      = D_N × w_sil,N × M_imp_met
N_core_old = w_core,N,old × M_core_old
N_core_new = N_core_old + N_met
```
である．

**Core δ¹⁵N：**
```
δ¹⁵N_core ← (δ¹⁵N_core × N_core_old + δ¹⁵N_met × N_met) / N_core_new
```

**Atmosphere δ¹⁵N：**
```
δ¹⁵N_atm ← δ¹⁵N_sil   （Δ¹⁵N^atm-sil = 0 より，Step 2 の平衡結果）
```
固化時の追加 outgas 項はなし（Step 3 で C_sil,N をマントルに保持するため）．

### Step 5：大気侵食（Section 6 参照）

ここで `m_a ≡ M_atm,total = ν_C M_atm,C + ν_H M_atm,H + ν_N M_atm,N` は atmosphere の総分子質量である．また `x_imp,i` は impactor 中の **元素 i の bulk mass fraction** とし，`x_imp,i = w_imp,i = C_imp,i × 10^-6` を用いる．

**GI フェーズ（GI1–GI9）：Schlichting (2015)（Section 6-B）：**
```
d        = M_imp / M_t                              （v_imp = v_esc より質量比のみ）
E_A_GI   = 0.4·d + 1.4·d² − 0.8·d³
ΔM_atm,i = −E_A_GI × M_atm,i
```

**暴走・寡占成長フェーズ（0.01→0.10）：Svetsov/Shuvalov 積分（Section 6-A）：**
```
ΔM_atm,i = −(E_p_ave × x_imp,i + E_A_ave × M_atm,i / m_a) × M_imp
```

δ¹⁵N_atm は侵食により変化しない（同位体中立；Section 6）．

### 部分溶融シナリオの自動処理

M_melt_sil は de Vries 式から衝突毎に計算され，全球溶融と部分溶融を統一的に扱う：

| シナリオ | M_melt_sil | 挙動 |
|---|---|---|
| 全球溶融（GI 大） | ≈ M_mantle | Sakuraba と同等 |
| 部分溶融（GI 小） | << M_mantle | melt pool のみ equilibration |

melt pool 外の固体 mantle は今回のステップでは平衡に参加しない．

---

## 6. 大気侵食モデル

Sakuraba のコード（Sakuraba.C）を参照し，フェーズごとに異なる侵食モデルを使用する．ζ, η は定数ではなく動的に計算される量であることに注意．

---

### 6-A. 暴走・寡占成長フェーズ（0.01→0.10）：Svetsov/Shuvalov 積分

Sakuraba.C の `f_etazeta_ave` に相当（Svetsov 2007; Shuvalov 2009 を速度・サイズ分布で数値積分）：

**離散化の規約**

- `0.01 → 0.10 M_Earth` は連続方程式をそのまま解かず，**固定質量刻みの precompute** として積分する
- デフォルト刻みは Sakuraba.C の `dA_accuracy = 10^-4` に合わせ，
  `dM_imp = 0.10 M_Earth × 10^-4 = 10^-5 M_Earth`
- したがって `0.01 → 0.10` は 9000 等間隔ステップで積分する
- 各ステップで現在の planet mass `M_t` を更新し，`R(M_t)`, `g(M_t)`, `v_esc(M_t)` を再計算したうえで Sections 3–5 を 1 回実行する
- 数値収束確認として，必要なら `dM_imp = 5×10^-6 M_Earth` に半減して最終 C/N/H と δ¹⁵N の差を比較する

```
dM_atm,i = −(E_p_ave × x_imp,i + E_A_ave × M_atm,i / m_a) × dM_imp
```

- **E_p_ave**（= ζ）：impactor 蒸気の惑星脱出効率（Shuvalov 2009 Eq.6 を積分）
- **E_A_ave**（= η）：大気直接侵食効率（Svetsov 2007 を積分）
- どちらも現ステップの大気スケール高度 H，密度 d_0，惑星半径 R_planet，脱出速度 v_esc の関数として**毎ステップ計算する**
- 積分範囲と分布は Sakuraba.C に合わせる
  - impactor size distribution：`dN/dD ∝ D^-2`（`f_mass2`），`D = 10^{3.5}–10^8 cm = 31.6 m – 1000 km`
  - impact velocity distribution：Rayleigh 型 `f(V) ∝ sqrt(V² − v_esc²) exp[−(V² − v_esc²)/(2v_esc²)]`，`V ∈ [v_esc, 4.324 v_esc]`
  - 数値積分は Simpson 法で行い，`E_A_ave`, `E_p_ave` の相対変化がともに `< 10^-2` で収束とする

---

### 6-B. GI フェーズ（GI1–GI9）：Schlichting (2015)

Sakuraba.C の `f_eta_GI` に相当（Sakuraba.C 1412–1416 行）：

```
d      = (V_imp × M_imp) / (v_esc × M_t)
E_A_GI = 0.4·d + 1.4·d² − 0.8·d³

ΔM_atm,i = −E_A_GI × M_atm,i   （現在の大気量の割合で侵食）
```

- `d` は衝突パラメータ（無次元）
- 我々のモデルでは V_imp = v_esc より **d = M_imp / M_t**（純粋な質量比）
- E_A_GI は既知量から計算される → 自由パラメータなし

---

### 6-C. 同位体中立性

衝突侵食は isotopically neutral：δ¹⁵N_atm は侵食により変化しない．

- 衝突侵食は tangent plane 上の大気を一括除去するバルクプロセス → ¹⁴N/¹⁵N を分別しない（Schlichting et al. 2015; Schlichting & Mukhopadhyay 2018）
- 仮に質量依存的逃脱を仮定しても，¹⁴N₂（28 amu）と ¹⁴N¹⁵N（29 amu）の質量差は ~3.6% にすぎず，δ¹⁵N への影響は数‰以下
- Shi et al. (2022)："atmospheric loss is bulk removal, isotopically neutral"
- 流体力学的散逸による ¹⁴N 優先逃脱は Xe 同位体制約から否定（Shi et al. 2022; Li et al. 2024）
