# egs (hg38): 0.85
# chip library size   5799920.0
# control library size   6907278.0
# bin size = 200
# gaps allowed = 3
# sults/0h/no_bigwig/0h_chip-W200-normalized.graph


# Find candidate islands exhibiting clustering ...
# /mnt/work/endrebak/software/anaconda/envs/py27/bin/python /home/endrebak/code/epic_paper/SICER/src/find_islands_in_pr.py -s hg38 -b /mnt/scratch/projects/epic_bencmarks/data/sicer_results/0h/no_bigwig/0h_chip-W200.graph -w 200 -g 600 -t 0.85 -e 1000 -f /mnt/scratch/projects/epic_bencmarks/data/sicer_results/0h/no_bigwig/0h_chip-W200-G600.scoreisland
# Species:  hg38
# Window_size:  200
# Gap size:  600
# E value is: 1000.0
# Total read count: 5799904.0
# Genome Length:  3088286401
# Effective genome Length:  2625043440
# Window average: 0.441890134969
# Window pvalue: 0.2
# ('self.average=', 0.44189013496858554)
# ('self.poisson_value[0]=', 0.6428202550436628)
# ('window_pvalue', 0.2)
# ('poisson', 0, 0.6428202550436628)
# ('poisson', 1, 0.28405592926178475)
# ('self.gap_contribution', 3.5822544757188606)
# ('self.boundary_contribution', 0.5447205833970497)
# Minimum num of tags in a qualified window:  2
# Generate the enriched probscore summary graph and filter the summary graph to get rid of ineligible windows
# Determine the score threshold from random background
# The score threshold is:  17.673
# Make and write islands
# Total number of islands:  38938


# Calculate significance of candidate islands using the control library ...
# /mnt/work/endrebak/software/anaconda/envs/py27/bin/python /home/endrebak/code/epic_paper/SICER/src/associate_tags_with_chip_and_control_w_fc_q.py -s hg38  -a /mnt/scratch/projects/epic_bencmarks/data/sicer_results/0h/no_bigwig/0h_chip-1-removed.bed -b /mnt/scratch/projects/epic_bencmarks/data/sicer_results/0h/no_bigwig/0h_input-1-removed.bed -d /mnt/scratch/projects/epic_bencmarks/data/sicer_results/0h/no_bigwig/0h_chip-W200-G600.scoreisland -f 150 -t 0.85 -o /mnt/scratch/projects/epic_bencmarks/data/sicer_results/0h/no_bigwig/0h_chip-W200-G600-islands-summary
# chip library size   5799920.0
# control library size   6907278.0
# Total number of chip reads on islands is:  1085652
# Total number of control reads on islands is:  383035


# Identify significant islands using FDR criterion ...
# /mnt/work/endrebak/software/anaconda/envs/py27/bin/python /home/endrebak/code/epic_paper/SICER/src/filter_islands_by_significance.py -i /mnt/scratch/projects/epic_bencmarks/data/sicer_results/0h/no_bigwig/0h_chip-W200-G600-islands-summary -p 1.0 -c 7 -o /mnt/scratch/projects/epic_bencmarks/data/sicer_results/0h/no_bigwig/0h_chip-W200-G600-islands-summary-FDR1.0
# Given significance 1.0 ,  there are 38938 significant islands


# Done!
import pytest

# from epic2.src.statistics import compute_background_probabilities
from epic2.src.SICER_stats import compute_score_threshold


def test_compute_background_probabilities():

    # 17.673
    result = compute_score_threshold(5799920, 200, 2625043440, 600)

    print(result)

    assert result[0] - 17.673 < 0.01



def test_compute_background_probabilities2():

    # 19.112
    result = compute_score_threshold(20227554, 200, 2625043440, 600)
    # result = compute_score_threshold(20227554, 200, 2625043440, 600)

    print(result)

    assert result[0] - 19.112 < 0.01


def test_compute_background_probabilities3():

    # 19.112
    result = compute_score_threshold(44590598, 200, 2625043440, 600)
    # result = compute_score_threshold(20227554, 200, 2625043440, 600)

    print(result)

    assert result[0] - 28.356 < 0.01
    assert result[1] == 6


def test_compute_background_probabilities4():

    # 19.112
    result = compute_score_threshold(44590598, 200, 2625043440, 600)
    # result = compute_score_threshold(20227554, 200, 2625043440, 600)

    print(result)

    assert result[0] - 28.356 < 0.01
    assert result[1] == 6
