#!/usr/bin/env python3
import subprocess
import sys
import os
import re

def run_cmd(cmd):
    print(f"正在执行: {cmd}")
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"命令执行失败:\n{res.stderr}", file=sys.stderr)
        sys.exit(1)
    return res.stdout

def load_fasta(fasta_path):
    sequences = {}
    name = None
    seq_lines = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    sequences[name] = "".join(seq_lines)
                name = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if name:
            sequences[name] = "".join(seq_lines)
    return sequences

def reverse_complement(seq):
    trans = str.maketrans("ATCGatcgNn", "TAGCtagcNn")
    return seq.translate(trans)[::-1]

def main():
    # 文件名配置
    ref_fasta = "CQ_chrY_v1.2_1gap.fasta"
    utg_fasta = "utg000737l.fa"
    output_fasta = "CQ_chrY_v1.2_filled.fasta"
    
    if not os.path.exists(ref_fasta) or not os.path.exists(utg_fasta):
        print("错误: 当前目录下未找到参考基因组或 utg 文件。", file=sys.stderr)
        sys.exit(1)

    paf_file = "alignment_utg.paf"
    # 使用 minimap2 进行比对，输出 PAF 格式
    run_cmd(f"minimap2 -x asm5 {ref_fasta} {utg_fasta} > {paf_file}")

    alns = []
    with open(paf_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            
            mapq = int(parts[11])
            # 过滤低质量比对，确保锚点唯一且可靠
            if mapq < 20: 
                continue
                
            aln = {
                'q_name': parts[0],
                'q_len': int(parts[1]),
                'q_start': int(parts[2]),
                'q_end': int(parts[3]),
                'strand': parts[4],
                't_name': parts[5],
                't_len': int(parts[6]),
                't_start': int(parts[7]),
                't_end': int(parts[8]),
                'matches': int(parts[9]),
                'aln_len': int(parts[10]),
                'mapq': mapq
            }
            alns.append(aln)

    if not alns:
        print("错误: 未找到 MapQ >= 20 的可靠比对锚点，无法执行填补。", file=sys.stderr)
        sys.exit(1)

    # 按照 (目标染色体, 链方向) 分组，寻找得分最高的主比对组
    from collections import defaultdict
    groups = defaultdict(list)
    for a in alns:
        groups[(a['t_name'], a['strand'])].append(a)
        
    # 取匹配碱基总数最多的一组作为填补依据
    best_group_key = max(groups.keys(), key=lambda k: sum(a['matches'] for a in groups[k]))
    best_alns = groups[best_group_key]
    
    # 按照参考基因组上的起始位置排序
    best_alns.sort(key=lambda a: a['t_start'])
    
    t_name, strand = best_group_key
    t_start_min = best_alns[0]['t_start']
    t_end_max = best_alns[-1]['t_end']
    
    # 无论正负链，PAF的 q_start 和 q_end 都是基于 forward 序列的坐标
    q_start_min = min(a['q_start'] for a in best_alns)
    q_end_max = max(a['q_end'] for a in best_alns)

    # 加载序列
    ref_seqs = load_fasta(ref_fasta)
    utg_seqs = load_fasta(utg_fasta)
    ref_seq = ref_seqs[t_name]
    utg_seq = utg_seqs[best_alns[0]['q_name']]

    # 统计 Gap (N) 信息
    replaced_ref_seq = ref_seq[t_start_min:t_end_max]
    ns = re.findall(r'[Nn]+', replaced_ref_seq)
    total_n_len = sum(len(n) for n in ns)

    # ================= 打印报告 =================
    print("\n==================================================")
    print(" [1] 替换位置与长度信息 (Replacement Summary)")
    print("==================================================")
    print(f"目标参考序列 (Target)    : {t_name}")
    print(f"参考序列替换区间         : {t_start_min} - {t_end_max}")
    print(f"被替换的参考总长度       : {t_end_max - t_start_min} bp")
    print(f"该区间包含的 Gap(N) 长度 : {total_n_len} bp")
    print("")
    print(f"插入序列 (Query Utg)     : {best_alns[0]['q_name']}")
    print(f"插入序列提取区间         : {q_start_min} - {q_end_max}")
    print(f"实际插入新序列长度       : {q_end_max - q_start_min} bp")
    print(f"比对方向 (Strand)        : {strand} （负链会自动反向互补替换）")

    print("\n==================================================")
    print(" [2] 上下游锚点序列与比对情况 (Anchors & Alignment)")
    print("==================================================")
    
    if len(best_alns) == 1:
        print("比对模式判断: Single-alignment (单片段直接跨越整个 Gap)")
        aln = best_alns[0]
        identity = (aln['matches'] / aln['aln_len']) * 100
        print(f"\n▶ 整体跨越锚点信息:")
        print(f"  - 锚定 Ref 区间: {aln['t_start']} - {aln['t_end']} (总跨度: {aln['t_end']-aln['t_start']} bp)")
        print(f"  - 匹配碱基数   : {aln['matches']} bp")
        print(f"  - 锚点比对质量 : MapQ = {aln['mapq']}")
        print(f"  - 整体比对一致率: {identity:.2f}%")
        
        # 尝试通过N的分布划分上下游锚点
        n_iters = list(re.finditer(r'[Nn]+', replaced_ref_seq))
        if n_iters:
            main_gap = max(n_iters, key=lambda m: m.end() - m.start())
            print(f"\n▶ 根据 N 序列分布推断的实际侧翼锚点长度:")
            print(f"  - 纯上游有效锚点长度: {main_gap.start()} bp")
            print(f"  - 纯下游有效锚点长度: {len(replaced_ref_seq) - main_gap.end()} bp")
            
    else:
        print(f"比对模式判断: Split-alignment (Gap 过大或序列差异导致比对断开为 {len(best_alns)} 个锚点片段)")
        up_aln = best_alns[0]
        down_aln = best_alns[-1]
        
        up_id = (up_aln['matches'] / up_aln['aln_len']) * 100
        down_id = (down_aln['matches'] / down_aln['aln_len']) * 100
        
        print(f"\n▶ 上游锚点 (Upstream Anchor):")
        print(f"  - Ref 区间     : {up_aln['t_start']} - {up_aln['t_end']} (长度: {up_aln['t_end']-up_aln['t_start']} bp)")
        print(f"  - 匹配碱基数   : {up_aln['matches']} bp")
        print(f"  - 锚点比对质量 : MapQ = {up_aln['mapq']}")
        print(f"  - 锚点比对一致率: {up_id:.2f}%")
        
        print(f"\n▶ 下游锚点 (Downstream Anchor):")
        print(f"  - Ref 区间     : {down_aln['t_start']} - {down_aln['t_end']} (长度: {down_aln['t_end']-down_aln['t_start']} bp)")
        print(f"  - 匹配碱基数   : {down_aln['matches']} bp")
        print(f"  - 锚点比对质量 : MapQ = {down_aln['mapq']}")
        print(f"  - 锚点比对一致率: {down_id:.2f}%")

    print("\n==================================================")
    print(" 开始执行序列替换...")
    
    # 提取替换序列
    replacement_seq = utg_seq[q_start_min:q_end_max]
    
    # 负链处理：无论 minimap2 报告片段如何，q_start_min/max 都是正链的物理坐标
    # 如果比对在负链，我们需要把这个区间内的 forward 序列取反向互补再放进去
    if strand == '-':
        replacement_seq = reverse_complement(replacement_seq)

    # 组装新序列
    new_chr_seq = ref_seq[:t_start_min] + replacement_seq + ref_seq[t_end_max:]
    ref_seqs[t_name] = new_chr_seq

    # 输出新文件
    with open(output_fasta, 'w') as out:
        for name, seq in ref_seqs.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

    print(f" 序列替换成功！\n 填补后的基因组已保存为: {output_fasta}")
    print("==================================================\n")

if __name__ == '__main__':
    main()