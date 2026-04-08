おすすめの順番はこれです。

まず struct を作る
PlanetState, ReservoirState, Impactor, StepInput, StepResult, Constants くらい。単位は最初に固定して、kg / ppm / K / GPa に統一します。

次に Section 2 のシナリオをコード化する
Impactor シナリオ (line 15) を std::vector<Impactor> に落として、0.01 -> 0.10 の precompute、GI1-9、late veneer を並べます。ここは物理計算なしでよくて、まず「どの順番で何が入るか」を固定します。

惑星の基本関数を作る
radius_from_mass(), gravity(), escape_velocity(), mantle_thickness()。これは後続の全部から使います。

melt pool 周りを実装する
3-2 (line 159), 3-3 (line 210) の V_melt -> d_pool -> P_eq -> T_eq を先に完成させます。ここは独立してテストしやすいです。

次に「簡易版の1衝突更新」を作る
Section 5 (line 648) の Mixing -> Equilibrium -> Reservoir update -> isotope -> erosion の流れだけ先に通します。
この段階では D_N, D_H, D_C, ΔIW_surface, 侵食効率は仮の定数でいいです。目的は「質量保存して最後まで1 run 回ること」です。

そのあと実データの分配係数を入れる
4-1 N (line 506), 4-2 H (line 560), 4-3 C (line 586) を差し込みます。D_N, Δ15N, D_H, D_C, X_O をここで完成させます。

self-oxidation を入れる
3-5 (line 296), 3-6 (line 414) は一番重いので後回しがいいです。
特に Fe3+/ΣFe, M_Fe0,oxid, ΔIW_surface は数値計算が多いので、まず本体が回る状態を作ってから追加した方がデバッグしやすいです。

最後に大気侵食を入れる
Section 6 (line 915) を最後に実装します。GI は比較的簡単ですが、0.01 -> 0.10 の Svetsov/Shuvalov 積分は重いので最後で十分です。

コードの呼び順は、最終的にはこう