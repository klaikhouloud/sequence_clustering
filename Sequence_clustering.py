# -*- coding: utf-8 -*-
import itertools
import logging
from Bio.SeqUtils.lcc import lcc_simp
from Bio.SeqUtils import gc_fraction  # Correction ici

def makelog(stri, do_print=True):
    if do_print:
        print(stri)
    logging.debug(stri)

def loadcluster(cluster_file):
    cluster_file, cluster_dic = open(cluster_file), {}
    cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
    for cluster in cluster_groups:
        name = next(cluster).strip()  # Correction Python 3
        seqs = [seq.split('>')[1].split('...')[0] for seq in next(cluster_groups)]
        cluster_dic[name[1:]] = seqs
    return cluster_dic

def filtercluster(cluster_dic, minimum, candidates):
    filtered_dic = {}
    for cluster in set(cluster_dic.keys()):
        if len(cluster_dic[cluster]) >= minimum:
            cluster_positions = {}
            for k in cluster_dic[cluster]:
                if isinstance(candidates[k]['start'], int) and isinstance(candidates[k]['end'], int):
                    if not candidates[k]['record'] in cluster_positions:
                        cluster_positions[candidates[k]['record']] = []
                    cluster_positions[candidates[k]['record']].append( (candidates[k]['start'], candidates[k]['end']) )
            total_len = 0
            for k,v in cluster_positions.items():
                overlapped_cluster = merge_overlap(v)
                total_len += len(overlapped_cluster)
            if total_len >= minimum:
                filtered_dic[cluster] = cluster_dic[cluster]
    return filtered_dic

def merge_overlap(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)
            else:
                merged.append(higher)
    return merged

def cluster2seq(cluster_dic, candidates, outfile):
    filter_file = open(outfile, 'w')
    last_cluster = None
    for cluster_id, cluster_seqs in cluster_dic.items():
        for seq_id in cluster_seqs:
            cand_id, description, sequence = candidates[seq_id]['id'], candidates[seq_id]['description'],candidates[seq_id]['seq']
            header = ">%s %s" % (cand_id, description)
            if not cluster_id == last_cluster:
                filter_file.write('{0}\n'.format('-'*10))
                last_cluster = cluster_id
            filter_file.write('{0}\n{1}\n'.format(header, '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])))
    filter_file.close()

def complex_enough(seq):
    complexity = lcc_simp(seq.upper())
    if complexity < 1.25:
        return False
    gc = gc_fraction(seq.upper()) * 100  # Correction ici : gc_fraction retourne une fraction
    if gc < 15 or gc > 95:
        return False
    return True

def cluster(file_names, candidates, min_copy_number, FSL, workers):
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
    from subprocess import Popen, PIPE
    from collections import OrderedDict
    import os, shutil
    import math

    makelog("Clustering")
    cmd_list = [
    './vsearch-2.7.1/bin/vsearch',
    '--cluster_fast',file_names['file_candidates_fasta'],
    #'--consout',file_names['file_representative'],
    '--threads',str(workers),
    '--strand','both',
    '--clusters',file_names['file_temp_cluster'],
    '--iddef','1',
    '--id', '0.8']
    makelog(' '.join(cmd_list))
    p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
    out,err = p.communicate()

    makelog("Clustering done")
    makelog("Filtering clusters")
    #count for minimum file length
    clusters_dic = {}
    list_dir = os.listdir(file_names['file_temp_cluster_dir'])
    makelog("Initial clusters: %i" % (len(list_dir),))
    for fn in list_dir:
        if os.path.isfile(file_names['file_temp_cluster_dir'] + fn):
            fh = open(file_names['file_temp_cluster_dir'] + fn)
            n = 0
            for line in fh:
                if line.startswith(">"):
                    n += 1
                    id_seq = line[1:line.find(" ")]
                    if fn in clusters_dic:
                        clusters_dic[fn].append(id_seq)
                    else:
                        clusters_dic[fn] = [id_seq]
            fh.close()
    
    filtered_clusters = filtercluster(clusters_dic, min_copy_number, candidates)
    unique_clusters = set(filtered_clusters.keys())
    num_clusters = len(unique_clusters)

    # Création de l'aligner une seule fois avant la boucle
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1

    for current_cluster in unique_clusters:
        candidates_in_cluster = filtered_clusters[current_cluster]
        new_min_copy_number = min_copy_number
        sum_diff_fs_cluster = 0
        for x in candidates_in_cluster:
            totally_different_fs = True
            cand_x = candidates[x]

            fs_right_1 = cand_x['fs_right']
            fs_left_1 = cand_x['fs_left']

            if fs_left_1 == '' or fs_right_1 == '' or not isinstance(fs_left_1,str) or not isinstance(fs_right_1,str):
                totally_different_fs = False
                continue
            if not complex_enough(fs_right_1) or not complex_enough(fs_left_1):
                totally_different_fs = False
                continue

            fs_right_1 = fs_right_1.upper()
            fs_left_1 = fs_left_1.upper()
            fs_right_1_plus_mite =  cand_x['seq'][-FSL:].upper() + fs_right_1
            fs_left_1_plus_mite = fs_left_1 +  cand_x['seq'][0:FSL].upper()

            at_least_one = False
            for y in candidates_in_cluster:
                cand_y = candidates[y]
                if cand_x['candidate_id'] == cand_y['candidate_id']:
                    continue
                fs_right_2 = cand_y['fs_right']
                fs_left_2 = cand_y['fs_left']

                if fs_right_2 == '' or fs_left_2 == '' or not isinstance(fs_right_2,str) or not isinstance(fs_left_2,str):
                    continue
                if not complex_enough(fs_right_2) or not complex_enough(fs_left_2):
                    continue

                fs_right_2 = fs_right_2.upper()
                fs_left_2 = fs_left_2.upper()
                fs_right_2_plus_mite =  cand_y['seq'][-FSL:].upper() + fs_right_2
                fs_left_2_plus_mite = fs_left_2 + cand_y['seq'][0:FSL].upper()

                fs_left_1_rc = Seq(fs_left_1).reverse_complement()
                fs_right_1_rc = Seq(fs_right_1).reverse_complement()

                fs_left_1_plus_mite_rc = Seq(fs_left_1_plus_mite).reverse_complement()
                fs_right_1_plus_mite_rc = Seq(fs_right_1_plus_mite).reverse_complement()

                # Utilisation de PairwiseAligner au lieu de pairwise2
                score_r1_r2 = aligner.score(fs_right_1, fs_right_2)
                score_l1_l2 = aligner.score(fs_left_1, fs_left_2)
                score_l1rc_r2 = aligner.score(fs_left_1_rc, fs_right_2)
                score_r1rc_l2 = aligner.score(fs_right_1_rc, fs_left_2)

                score_r1_r2_plus_mite = aligner.score(fs_right_1, fs_right_2_plus_mite)
                score_l1_l2_plus_mite = aligner.score(fs_left_1, fs_left_2_plus_mite)
                score_l1rc_r2_plus_mite = aligner.score(fs_left_1_rc, fs_right_2_plus_mite)
                score_r1rc_l2_plus_mite = aligner.score(fs_right_1_rc, fs_left_2_plus_mite)

                max_score = max(score_r1_r2,score_l1_l2,score_l1rc_r2,score_r1rc_l2,
                                score_r1_r2_plus_mite,score_l1_l2_plus_mite,
                                score_l1rc_r2_plus_mite,score_r1rc_l2_plus_mite)

                if max_score == []:
                    max_score = 0
                max_score /= FSL
                at_least_one = True
                if max_score > 0.5:
                    totally_different_fs = False
                    break
               
            if totally_different_fs and at_least_one:
                sum_diff_fs_cluster += 1
            if sum_diff_fs_cluster >= new_min_copy_number:
                break

        if sum_diff_fs_cluster < new_min_copy_number:
            del filtered_clusters[current_cluster]

    ordered_cluster = OrderedDict(sorted(filtered_clusters.items(), key=lambda t: t[0]))

    makelog("Clusters: " + str(len(filtered_clusters)))
    buffer_rec = []
    count = 1
    family_number = 1
    buffer_nr = []
    print("Clés disponibles dans file_names :", list(file_names.keys()))
    output_file = file_names['file_family_fasta']
    with open(output_file, 'w') as output_handle:
        for cluster_id, members in ordered_cluster.items():
            for member in members:
                seq = candidates[member]['seq']
                rec = SeqRecord(Seq(seq), id=str(count), description=str(family_number))
                buffer_rec.append(rec)
                count += 1
            family_number += 1
        SeqIO.write(buffer_rec, output_handle, 'fasta')

    return filtered_clusters

