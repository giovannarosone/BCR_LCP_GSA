The output format for eBWT/LCP/GSA(DA/SA) can be:
- EGSA (OUTPUT_FORMAT must be 1). The end-marker in .bwt file is the symbol '\0' [https://github.com/felipelouza/egsa]
- at most four files: .ebwt, .lcp, .da, posSA (OUTPUT_FORMAT != 1). The end-marker in .bwt file is the symbol '#'
    - if OUTPUT_FORMAT == 0, the output format of BCR is at most 4 files - built one after the other
    - if OUTPUT_FORMAT == 1, the output format of BCR is as the output of EGSA (.gesa file). BUILD_LCP, BUILD_DA and BUILD_SA must be set to 1. Please, set the types as in eGSA
    - if OUTPUT_FORMAT == 2, the output format of BCR is a unique file .egsa. BUILD_LCP must be set to 1 (we do not use a struct), BUILD_DA and BUILD_SA could be set to either 0 or 1.  Order: ebwt, lcp, da, sa
    - if OUTPUT_FORMAT == 3, the output format of BCR is at most 4 files at the same time
    - if OUTPUT_FORMAT == 4, the output format of BCR is at most 3 files (ebwt, da), lcp, sa
    - if OUTPUT_FORMAT == 5, the output format of BCR is at most 3 files ebwt, (lcp, da), sa
    - if OUTPUT_FORMAT == 6, the output format of BCR is at most 3 files ebwt, lcp, (sa, da)
