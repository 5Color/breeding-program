#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CLI 품종개량 시뮬레이터 v7 + 생물학 제약 (코돈/GC/stop/AA/마커)

- 초기 염기서열을 20bp 이상이면 어떤 길이든 그대로 사용 (상한 없음)
- 초기 DNA 미제공 시, 사용자가 원하는 DNA 길이(≥20)를 입력받아 랜덤 시작
- 단위 규칙: Sweetness=Brix(%), 숙성/저장/조기재배=days, 그 외 score/10
- 인덱스로 작물 선택, 형질 목표값 직접 입력, tqdm 진행바, 최종 염기서열 출력
- ✅ 초기 염기서열(init)이 주어지면 다음 “생물학 제약”이 추가로 적용됨:
    * 코돈(3염기) 단위로만 생물학적 평가
    * 원본 대비 아미노산 서열이 많이 바뀌면 패널티
    * 새로 생긴 stop codon(TAA/TAG/TGA)에 강한 패널티
    * 원본 GC%와 너무 다르면 패널티
    * 마커 코돈 위치에 ‘고GC 동의 코돈’이면 보너스
"""

import sys, random, argparse, textwrap
from typing import List, Tuple, Dict, Any

from tqdm import tqdm


# -----------------------------
# 작물/형질 정의
# (name, weight, goal, base, icon, unit, vmin, vmax)
# -----------------------------
CROPS: Dict[str, List[Tuple[Any, ...]]] = {
    "Tomato": [
        ("수확량(Yield)",           0.25, 7.0,  5.0, "🍅", "score/10",     0.0, 10.0),
        ("내병성(Resistance)",      0.20, 8.0,  5.0, "🌿", "score/10",     0.0, 10.0),
        ("당도(Sweetness)",         0.15,12.0,  8.0, "🍓", "Brix(%)",      4.0, 25.0),
        ("숙성속도(Maturing)",      0.15,15.0, 30.0, "⏳", "days↓(fast)",  5.0, 60.0),
        ("영양가(Nutrition)",       0.15, 7.0,  5.0, "💪", "score/10",     0.0, 10.0),
        ("가뭄내성(DroughtTol.)",   0.10, 7.0,  5.0, "🏜️", "score/10",     0.0, 10.0),
    ],
    "Pepper": [
        ("매운맛(Spiciness)",       0.25, 8.0,  5.0, "🌶️", "score/10",     0.0, 10.0),
        ("크기(Size)",              0.20, 6.0,  5.0, "📏", "score/10",      0.0, 10.0),
        ("색상선명도(Color)",       0.15, 7.0,  5.0, "🎨", "score/10",      0.0, 10.0),
        ("내병성(Resistance)",      0.15, 8.0,  5.0, "🛡️", "score/10",     0.0, 10.0),
        ("저장성(ShelfLife)",       0.15,25.0, 15.0, "📦", "days",          5.0, 60.0),
        ("수확량(Yield)",           0.10, 7.0,  5.0, "🫑", "score/10",     0.0, 10.0),
    ],
    "Rice": [
        ("수확량(Yield)",           0.25, 7.0,  5.0, "🌾", "score/10",     0.0, 10.0),
        ("병충해저항성(PestRes.)",  0.20, 8.0,  5.0, "🐛", "score/10",     0.0, 10.0),
        ("단백질(Protein)",         0.15, 6.0,  5.0, "🍚", "score/10",      0.0, 10.0),
        ("조기재배(EarlyCult.)",    0.15,20.0, 45.0, "☀️", "days↓(fast)", 10.0, 60.0),
        ("가뭄내성(DroughtTol.)",   0.15, 7.0,  5.0, "💧", "score/10",     0.0, 10.0),
        ("미질(GrainQuality)",      0.10, 8.0,  5.0, "✨", "score/10",      0.0, 10.0),
    ],
}
CROP_NAMES = list(CROPS.keys())
NUM_TRAITS = 6
BASES = ("A","C","G","T")

# -----------------------------
# 공통 유틸
# -----------------------------
def clean_dna(s: str) -> str:
    s = (s or "").upper()
    return "".join(ch for ch in s if ch in "ACGT")

def ask(prompt: str, cast, default=None, cond=lambda x: True):
    while True:
        raw = input(f"{prompt}{' ['+str(default)+']' if default is not None else ''}: ").strip()
        if raw == "" and default is not None:
            val = default
        else:
            try:
                val = cast(raw)
            except Exception:
                print("잘못된 입력입니다. 다시 입력해 주세요.")
                continue
        if cond(val):
            return val
        print("값 범위를 확인하세요.")

def ask_index(prompt: str, options: List[str], default_idx: int = 1) -> int:
    while True:
        print(prompt)
        for i, name in enumerate(options, 1):
            print(f"  {i}. {name}")
        raw = input(f"번호를 입력하세요 [{default_idx}]: ").strip()
        if raw == "": raw = str(default_idx)
        if raw.isdigit():
            idx = int(raw)
            if 1 <= idx <= len(options):
                return idx - 1
        print("잘못된 번호입니다. 다시 입력해 주세요.")

def even_split_ranges(n: int, k: int) -> List[Tuple[int,int]]:
    base = n // k
    r = n % k
    seglens = [(base + 1 if i < r else base) for i in range(k)]
    ranges = []
    cur = 0
    for L in seglens:
        ranges.append((cur, cur + L))
        cur += L
    return ranges

def rand_dna(n: int) -> str:
    return "".join(random.choice(BASES) for _ in range(n))

# -----------------------------
# 🔬 생물학 제약용 유틸 (코돈/GC/AA/마커)
# -----------------------------
def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return gc / len(seq)

def chunk_codons(seq: str) -> List[str]:
    n = (len(seq) // 3) * 3
    return [seq[i:i+3] for i in range(0, n, 3)]

def codons_to_seq(codons: List[str]) -> str:
    return "".join(codons)

# 표준 유전 암호
CODON_TABLE: Dict[str, str] = {
    # Phe
    "TTT": "F", "TTC": "F",
    # Leu
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Ile
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Met (Start)
    "ATG": "M",
    # Val
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Ser
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    # Pro
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Thr
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Ala
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyr
    "TAT": "Y", "TAC": "Y",
    # His
    "CAT": "H", "CAC": "H",
    # Gln
    "CAA": "Q", "CAG": "Q",
    # Asn
    "AAT": "N", "AAC": "N",
    # Lys
    "AAA": "K", "AAG": "K",
    # Asp
    "GAT": "D", "GAC": "D",
    # Glu
    "GAA": "E", "GAG": "E",
    # Cys
    "TGT": "C", "TGC": "C",
    # Trp
    "TGG": "W",
    # Arg
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGA": "R", "AGG": "R",
    # Ser (AG-)
    "AGT": "S", "AGC": "S",
    # Gly
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    # Stop
    "TAA": "*", "TAG": "*", "TGA": "*",
}
STOP_CODONS = {"TAA", "TAG", "TGA"}

AA_TO_CODONS: Dict[str, List[str]] = {}
for codon, aa in CODON_TABLE.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)

def translate_codons(codons: List[str]) -> str:
    aas = []
    for c in codons:
        aas.append(CODON_TABLE.get(c, "X"))
    return "".join(aas)

def pick_markers(num_codons: int, n_markers: int = 15) -> List[int]:
    if n_markers > num_codons:
        n_markers = num_codons
    return sorted(random.sample(range(num_codons), n_markers))

def choose_high_gc_codons_for_aa(aa: str) -> List[str]:
    codons = AA_TO_CODONS.get(aa, [])
    if not codons:
        return []
    scored = sorted(codons, key=lambda c: gc_content(c), reverse=True)
    best_gc = gc_content(scored[0])
    best = [c for c in scored if gc_content(c) == best_gc]
    return best

def build_favorable_codons(
    ref_codons: List[str], markers: List[int]
) -> Dict[int, List[str]]:
    favorable: Dict[int, List[str]] = {}
    ref_aa = translate_codons(ref_codons)
    for idx in markers:
        if idx >= len(ref_codons):
            continue
        aa = ref_aa[idx]
        if aa in ("*", "X"):
            continue
        best_codons = choose_high_gc_codons_for_aa(aa)
        if best_codons:
            favorable[idx] = best_codons
    return favorable

def bio_score(dna: str, bio: Dict[str, Any]) -> float:
    """
    생물학 제약 기반 추가 점수 (양수=보너스, 음수=패널티).
    """
    ref_codons: List[str] = bio["ref_codons"]
    ref_aa: str = bio["ref_aa"]
    ref_gc: float = bio["ref_gc"]
    favorable_codons: Dict[int, List[str]] = bio["favorable_codons"]

    stop_penalty: float = bio.get("stop_penalty", 40.0)
    aa_change_penalty: float = bio.get("aa_change_penalty", 30.0)
    gc_penalty_scale: float = bio.get("gc_penalty_scale", 100.0)
    marker_bonus: float = bio.get("marker_bonus", 5.0)

    score = 0.0

    codons = chunk_codons(dna)
    if not codons or not ref_codons:
        return 0.0
    seq_aa = translate_codons(codons)

    # 1) 마커 유리 코돈 보너스
    for idx, fav_list in favorable_codons.items():
        if idx < len(codons) and codons[idx] in fav_list:
            score += marker_bonus

    # 2) 아미노산 변경 비율 패널티
    min_len = min(len(seq_aa), len(ref_aa))
    if min_len > 0:
        diff = sum(1 for i in range(min_len) if seq_aa[i] != ref_aa[i])
        aa_diff_ratio = diff / min_len
        score -= aa_diff_ratio * aa_change_penalty

    # 3) stop codon 새로 생성 패널티
    max_len_cod = min(len(codons), len(ref_codons))
    for i in range(max_len_cod):
        c = codons[i]
        r = ref_codons[i]
        if c in STOP_CODONS and r not in STOP_CODONS:
            score -= stop_penalty

    # 4) GC% 유지 패널티
    cur_gc = gc_content(dna)
    diff_gc = abs(cur_gc - ref_gc)
    if diff_gc > 0.05:  # ±5%까지는 허용
        score -= (diff_gc - 0.05) * gc_penalty_scale

    return score

# -----------------------------
# 형질 계산 & 적합도
# -----------------------------
def seg_genotype_score(seg: str) -> float:
    if not seg:
        return 0.0
    m = {"A":0,"C":1,"G":2,"T":3}
    s = sum(m[ch] for ch in seg)
    return s / float(3 * len(seg))  # 0~1

def trait_value(traits, env_sev: float, seg: str, trait_idx: int):
    name, weight, goal, base_pheno, _icon, unit, vmin, vmax = traits[trait_idx]
    geno01 = seg_genotype_score(seg)      # 0~1
    genotype_score = geno01 * 100.0       # 0~100
    max_contrib = max(goal, base_pheno) + (vmax - vmin) * 0.5
    phenotype = base_pheno + geno01 * (max_contrib - base_pheno)
    noise = (random.random() - 0.5) * env_sev * ((vmax - vmin) * 0.05)
    phenotype = max(vmin, min(vmax, phenotype + noise))
    return genotype_score, phenotype

def fitness_of(
    dna: str,
    traits,
    env_sev: float,
    seg_ranges: List[Tuple[int,int]],
    bio: Dict[str, Any] = None
):
    """
    total_fitness = 형질 적합도(0~100 근처) + 생물학 점수(보너스/패널티)
    """
    total = 0.0
    details = []
    for i in range(NUM_TRAITS):
        a, b = seg_ranges[i]
        seg = dna[a:b]
        gscore, pheno = trait_value(traits, env_sev, seg, i)
        name, weight, goal, _base, icon, unit, vmin, vmax = traits[i]
        span = max(1e-6, (vmax - vmin))
        diff_std = abs(pheno - goal) / span
        trait_fit = max(0.0, 100.0 - 100.0 * diff_std)
        total += trait_fit * weight
        details.append({
            "trait": name, "icon": icon, "unit": unit, "goal": goal,
            "genotype": gscore, "phenotype": pheno, "fitness": trait_fit, "seg": seg
        })

    # 🔬 생물학 제약 점수 추가
    if bio is not None:
        total += bio_score(dna, bio)

    return total, details

# -----------------------------
# GA 연산
# -----------------------------
def roulette(pop):
    total = sum(p["fitness"] for p in pop)
    if total <= 0:
        return random.choice(pop)
    r = random.random() * total
    s = 0.0
    for ind in pop:
        s += ind["fitness"]
        if s >= r: return ind
    return pop[-1]

def crossover(dna1: str, dna2: str) -> str:
    if len(dna1) != len(dna2): raise ValueError("DNA length mismatch")
    if len(dna1) < 2: return dna1
    x = random.randint(1, len(dna1)-1)
    return dna1[:x] + dna2[x:]

def mutate(dna: str, rate: float=0.01) -> str:
    arr = list(dna)
    for i in range(len(arr)):
        if random.random() < rate:
            arr[i] = random.choice(("A","C","G","T"))
    return "".join(arr)

def evolve(pop_size: int, generations: int, traits, env_sev: float,
           dna_len: int, seg_ranges: List[Tuple[int,int]], init_dna: str=None,
           mutation_rate: float=0.01, preview_every: int=10,
           bio: Dict[str, Any] = None):
    population = []
    if init_dna:
        f, d = fitness_of(init_dna, traits, env_sev, seg_ranges, bio=bio)
        population.append({"dna": init_dna, "fitness": f, "details": d})
    while len(population) < pop_size:
        dna = rand_dna(dna_len)
        f, d = fitness_of(dna, traits, env_sev, seg_ranges, bio=bio)
        population.append({"dna": dna, "fitness": f, "details": d})

    best = max(population, key=lambda x: x["fitness"])
    bar = tqdm(range(1, generations+1), desc="세대 진화", ncols=80)
    for gen in bar:
        new_pop = []
        elite = max(population, key=lambda x: x["fitness"])
        new_pop.append(elite)

        while len(new_pop) < pop_size:
            p1 = roulette(population)
            p2 = roulette(population)
            child_dna = crossover(p1["dna"], p2["dna"])
            child_dna = mutate(child_dna, mutation_rate)
            f, d = fitness_of(child_dna, traits, env_sev, seg_ranges, bio=bio)
            new_pop.append({"dna": child_dna, "fitness": f, "details": d})

        population = new_pop
        cand = max(population, key=lambda x: x["fitness"])
        if cand["fitness"] > best["fitness"]:
            best = cand

        bar.set_postfix_str(f"best={best['fitness']:.2f}")
        if preview_every and gen % preview_every == 0:
            show = best['dna'][:60] + ("..." if len(best['dna'])>60 else "")
            bar.write(f"[{gen:>4}세대] best_fit={best['fitness']:.2f} dna={show}")

    return best

# -----------------------------
# 출력
# -----------------------------
def print_result(best, crop_name, seg_ranges: List[Tuple[int,int]], bio: Dict[str, Any] = None):
    print("\n" + "="*78)
    print(f"✅ 최적 개체 결과 ({crop_name})")
    print("-"*78)
    print(f"최적 적합도(형질+생물학): {best['fitness']:.2f}")
    print(f"최종 최적 염기서열({len(best['dna'])}bp):")
    seq = best['dna']
    for i in range(0, len(seq), 80):
        print(seq[i:i+80])
    print("-"*78)
    print("형질별 현상/목표 및 세그먼트")
    print("-"*78)
    print(f"{'Trait':22} {'Goal(unit)':>16} {'Pheno':>12} {'Fit':>7} {'Seg(len)':>14}  Genotype(0~100)")
    for i, t in enumerate(best["details"]):
        a,b = seg_ranges[i]
        seg = t["seg"]
        seg_label = f"{seg}({b-a})" if len(seg) <= 14 else f"{seg[:11]}...({b-a})"
        goal_unit = f"{t['goal']:.2f} {t['unit']}"
        print(f"{t['icon']} {t['trait'][:18]:18} "
              f"{goal_unit:>16} {t['phenotype']:12.2f} {t['fitness']:7.2f} "
              f"{seg_label:>14}  {t['genotype']:6.1f}")
    print("="*78)

    # 🔬 생물학 요약 추가 출력
    if bio is not None:
        print("\n[생물학 제약 요약]")
        ref_gc = bio["ref_gc"]
        best_gc = gc_content(best["dna"])
        print(f"원본 GC%: {ref_gc*100:.2f}% | 최종 GC%: {best_gc*100:.2f}%")

        ref_codons = bio["ref_codons"]
        best_codons = chunk_codons(best["dna"])
        ref_aa = bio["ref_aa"]
        best_aa = translate_codons(best_codons)

        min_len = min(len(ref_aa), len(best_aa))
        diff = sum(1 for i in range(min_len) if ref_aa[i] != best_aa[i])
        aa_diff_ratio = (diff / min_len) if min_len > 0 else 1.0
        print(f"아미노산 변경 비율: {aa_diff_ratio*100:.2f}% (len={min_len})")

        favorable_codons = bio["favorable_codons"]
        marker_total = len(favorable_codons)
        hit = 0
        for idx, fav_list in favorable_codons.items():
            if idx < len(best_codons) and best_codons[idx] in fav_list:
                hit += 1
        print(f"마커 유리 코돈 달성: {hit}/{marker_total}")

        # 일부 마커 상세
        print("\n[마커 상세 (최대 10개)]")
        shown = 0
        for idx, fav_list in favorable_codons.items():
            if idx >= len(ref_codons) or idx >= len(best_codons):
                continue
            ref_c = ref_codons[idx]
            best_c = best_codons[idx]
            ok = "OK" if best_c in fav_list else "--"
            print(f" idx {idx:4d}: ref={ref_c} -> best={best_c} [{ok}]")
            shown += 1
            if shown >= 10:
                break
        print("="*78)

    print("\n")

# -----------------------------
# 목표값 입력 (인터랙티브)
# -----------------------------
def prompt_goals(traits: List[Tuple[Any, ...]]) -> List[Tuple[Any, ...]]:
    new_traits = []
    print("\n--- 형질 목표값 설정 (빈칸=기본 유지) ---")
    for name, weight, goal, base, icon, unit, vmin, vmax in traits:
        while True:
            s = input(f"{icon} {name} 목표 ({unit}, 권장 {vmin}~{vmax}, 기본 {goal}): ").strip()
            if s == "":
                new_goal = goal
                break
            try:
                val = float(s)
                if val < vmin or val > vmax:
                    print(f"범위를 벗어났습니다. {vmin}~{vmax} 사이로 입력하세요.")
                    continue
                new_goal = val
                break
            except:
                print("숫자로 입력하세요.")
        new_traits.append((name, weight, new_goal, base, icon, unit, vmin, vmax))
    return new_traits

# -----------------------------
# 메인
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="품종개량 CLI v7 + 생물학 제약 (초기 DNA≥20bp면 그대로 사용, 상한 없음)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        예)
          python ga_breeding_cli_v7_bio.py
          python ga_breeding_cli_v7_bio.py --noninteractive --crop_idx 1 --pop 600 --gens 400 --env 2.0 --dna_len 120
          python ga_breeding_cli_v7_bio.py --noninteractive --crop_idx 2 --init ATGCA...TGAC   # 초기염기서열 기준 생물학 제약
        """),
    )
    parser.add_argument("--noninteractive", action="store_true")
    parser.add_argument("--crop_idx", type=int, default=1, help=f"작물 인덱스(1~{len(CROP_NAMES)})")
    parser.add_argument("--pop", type=int, default=500)
    parser.add_argument("--gens", type=int, default=400)
    parser.add_argument("--env", type=float, default=2.0)
    parser.add_argument("--init", type=str, default="")
    parser.add_argument("--dna_len", type=int, default=24)
    parser.add_argument("--mut", type=float, default=0.01)
    parser.add_argument("--preview", type=int, default=10)
    args = parser.parse_args()

    if not (1 <= args.crop_idx <= len(CROP_NAMES)):
        print(f"--crop_idx는 1~{len(CROP_NAMES)} 사이여야 합니다.", file=sys.stderr); sys.exit(1)

    crop_name = CROP_NAMES[args.crop_idx - 1]
    base_traits = CROPS[crop_name]

    if not args.noninteractive:
        print("=== 품종개량 CLI v7 + 생물학 제약 ===")
        idx = ask_index("작물을 선택하세요:", CROP_NAMES, default_idx=args.crop_idx)
        crop_name = CROP_NAMES[idx]
        base_traits = CROPS[crop_name]

        pop  = ask("개체군 크기(pop_size)", int, default=args.pop,  cond=lambda x: x>=10)
        gens = ask("세대 수(generations)", int, default=args.gens, cond=lambda x: x>=1)
        env  = ask("환경강도(env_severity, 0~5)", float, default=args.env, cond=lambda x: 0<=x<=5)

        traits = prompt_goals(base_traits)

        # ✅ 초기 염기서열: 20bp 이상이면 어떤 길이든 그대로 사용
        init_raw = input("\n초기 염기서열(선택, 길이≥20, A/C/G/T): ").strip().upper()
        init = clean_dna(init_raw)
        if init and len(init) < 20:
            print("❌ 초기 염기서열은 최소 20bp여야 합니다. 무시하고 길이 입력 단계로 넘어갑니다.")
            init = ""
        if init:
            dna_len = len(init)  # 상한 없음, 그대로 진행
            print(f"초기 염기서열 길이 감지: {dna_len}bp (이 길이에 맞춰 진행)")
        else:
            dna_len = ask("DNA 길이(랜덤 시작, 길이≥20)", int,
                          default=max(24, args.dna_len), cond=lambda x: x>=20)

        mut  = ask("돌연변이율(염기당)", float, default=args.mut, cond=lambda x: 0<=x<=0.5)
        preview = ask("미리보기 간격(세대)", int, default=args.preview, cond=lambda x: x>=0)
    else:
        pop, gens, env, mut, preview = args.pop, args.gens, args.env, args.mut, args.preview
        traits = base_traits
        init = clean_dna(args.init)
        if init and len(init) < 20:
            print("[경고] 초기 염기서열은 최소 20bp입니다. 무시하고 랜덤 시작.", file=sys.stderr)
            init = ""
        dna_len = len(init) if init else max(20, args.dna_len)

    seg_ranges = even_split_ranges(dna_len, NUM_TRAITS)

    # 🔬 생물학 제약 설정 (초기 염기서열이 있을 때만)
    bio = None
    if init:
        ref_dna = init
        ref_gc = gc_content(ref_dna)
        ref_codons = chunk_codons(ref_dna)
        ref_aa = translate_codons(ref_codons)
        if len(ref_codons) > 0:
            markers = pick_markers(len(ref_codons), n_markers=15)
            favorable_codons = build_favorable_codons(ref_codons, markers)
            bio = {
                "ref_gc": ref_gc,
                "ref_codons": ref_codons,
                "ref_aa": ref_aa,
                "favorable_codons": favorable_codons,
                # 하이퍼파라미터 (원하면 여기 값 조정해서 실험 가능)
                "stop_penalty": 40.0,
                "aa_change_penalty": 30.0,
                "gc_penalty_scale": 100.0,
                "marker_bonus": 5.0,
            }
            print("\n[생물학 제약 활성화]")
            print(f"- 원본 길이: {len(ref_dna)}bp, 코돈 수: {len(ref_codons)}")
            print(f"- 원본 GC%: {ref_gc*100:.2f}%")
            print(f"- 마커 코돈 개수: {len(favorable_codons)} (stop/알수없는 위치 제외)")
        else:
            print("\n[주의] 초기 염기서열이 있지만 코돈이 없습니다. (길이<3?) 생물학 제약 비활성.")

    random.seed()
    best = evolve(pop, gens, traits, env, dna_len, seg_ranges,
                  init_dna=init, mutation_rate=mut, preview_every=preview,
                  bio=bio)
    print_result(best, crop_name, seg_ranges, bio=bio)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[중단됨]")
