# use rule fa_to_bam as chrr_fa_to_bam with:
#     input:
#         fa = config['ref']['fa_chrr']
#     output:
#         bam = config['ref']['fa_chrr_bam']
#
# use rule index_bam as chrr_ind_bam with:
#     input:
#         bam = rules.chrr_fa_to_bam.output.bam
#     output:
#         ind = config['ref']['fa_chrr_bam_ind']

# use rule compute_gc as chrr_gc_bw with:
#     input:
#         fa = config['ref']['fa_chrr'],
#         fa_ind = config['ref']['fa_chrr_ind']
#     output:
#         bw = config['ref']['fa_chrr_gc_bw']
