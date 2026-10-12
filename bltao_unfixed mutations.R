rm(list = ls())

# ============================================================
# 0) 加载包
# ============================================================

library(tidyr)
library(Biostrings)
library(readxl)
library(dplyr)
library(stringr)
library(janitor)
library(readr)
library(purrr)
library(tibble)
library(openxlsx)
library(ggplot2)

# ============================================================
# 1) 路径和文件名
# ============================================================

setwd("/home/data/t040512/UnfixedMutation")

mut_xlsx <- "high_confidence_loci_simplified.xlsx"
gb_file  <- "H37Rv.gb"
gene_sites_rds <- "H37Rv_gene_sites.rds"

excel_file <- "unfixed_mutation_analysis_results.xlsx"

# ============================================================
# 2) 读取变异表，并正确判断 synonymous / nonsynonymous
# ============================================================

mut <- read_excel(mut_xlsx) |>
  clean_names() |>
  transmute(
    sample       = as.character(sample),
    position     = as.numeric(position),
    refvar       = str_trim(as.character(ref_var)),
    freq         = suppressWarnings(
      as.numeric(str_replace_all(as.character(frequency), "%", ""))
    ),
    depth        = as.numeric(depth),
    mutationtype = str_to_lower(str_trim(as.character(mutation_type))),
    gene_id      = str_extract(str_trim(as.character(gene)), "Rv\\d{4}[A-Za-z]?")
  ) |>
  mutate(
    freq = ifelse(freq > 1, freq / 100, freq),
    
    mutationtype_norm = mutationtype |>
      str_to_lower() |>
      str_trim() |>
      str_replace_all("[-_ ]", ""),
    
    is_nonsyn = mutationtype_norm %in% c(
      "nonsynonymous",
      "nonsynonymousvariant",
      "nonsynonymousmutation",
      "missense",
      "missensevariant",
      "missensemutation",
      "n"
    ),
    
    is_syn = mutationtype_norm %in% c(
      "synonymous",
      "synonymousvariant",
      "synonymousmutation",
      "silent",
      "silentvariant",
      "silentmutation",
      "s"
    )
  ) |>
  mutate(
    # ============================================================
    # 新增：Intergenic / intergenic 非固定性计数口径
    # ============================================================
    is_intergenic = mutationtype_norm %in% c(
      "intergenic",
      "intergenicvariant",
      "intergenicmutation",
      "ig",
      "intergenic_variant",
      "intergenic_mutation"
    )
  ) |>
  filter(!is.na(sample), !is.na(freq), !is.na(depth)) |>
  filter(depth >= 20, freq >= 0.01, freq <= 0.90)

# ============================================================
# 3) 检查 mutation_type 分类
# ============================================================

cat("\n---- mutation_type 原始分布 ----\n")
print(table(mut$mutationtype, useNA = "ifany"))

cat("\n---- mutation_type 判断结果 ----\n")
print(
  mut |>
    group_by(mutationtype, mutationtype_norm, is_nonsyn, is_syn, is_intergenic) |>
    summarise(n = n(), .groups = "drop") |>
    arrange(desc(n))
)

cat(sprintf(
  "\nN(nonsyn) = %d\nS(syn)    = %d\nI(inter)  = %d\n",
  sum(mut$is_nonsyn, na.rm = TRUE),
  sum(mut$is_syn,    na.rm = TRUE),
  sum(mut$is_intergenic, na.rm = TRUE)
))

if (any(mut$is_nonsyn & mut$is_syn, na.rm = TRUE)) {
  stop("错误：存在同时被判定为 nonsynonymous 和 synonymous 的突变。")
}

# （可选但保留你原逻辑风格：不检查 intergenic 是否与 nonsyn/syn 重叠）

cat("\n---- synonymous 检查 ----\n")
print(
  mut |>
    filter(mutationtype_norm == "synonymous") |>
    group_by(is_nonsyn, is_syn, is_intergenic) |>
    summarise(n = n(), .groups = "drop") |>
    arrange(desc(n))
)

cat("\n---- nonsynonymous 检查 ----\n")
print(
  mut |>
    filter(mutationtype_norm == "nonsynonymous") |>
    group_by(is_nonsyn, is_syn, is_intergenic) |>
    summarise(n = n(), .groups = "drop") |>
    arrange(desc(n))
)

cat("\n---- 未识别为 synonymous / nonsynonymous 的 mutation_type ----\n")
print(
  mut |>
    filter(!is_syn, !is_nonsyn) |>
    group_by(mutationtype, mutationtype_norm) |>
    summarise(n = n(), .groups = "drop") |>
    arrange(desc(n))
)

# ============================================================
# 新增：每个样本的非固定性突变数量汇总（Nonsyn/Syn/Intergenic）
# 口径：非固定性 = Nonsynonymous + Synonymous + Intergenic
# ============================================================

per_sample_nonfixed_variant_counts <- mut |>
  summarise(
    nonsynonymous_nonfixed = sum(is_nonsyn, na.rm = TRUE),
    synonymous_nonfixed    = sum(is_syn, na.rm = TRUE),
    intergenic_nonfixed    = sum(is_intergenic, na.rm = TRUE),
    total_nonfixed        = nonsynonymous_nonfixed + synonymous_nonfixed + intergenic_nonfixed,
    .by = sample
  ) |>
  arrange(desc(total_nonfixed))

write_tsv(per_sample_nonfixed_variant_counts,
          "per_sample_nonfixed_variant_counts.tsv")

cat("\n---- 每个样本 非固定性突变数量（前 10 行）----\n")
print(head(per_sample_nonfixed_variant_counts, 10))

# ============================================================
# 4) 解析 GenBank 并计算每个基因的 N_sites / S_sites
# ============================================================

if (file.exists(gene_sites_rds)) {
  
  cat(sprintf("\n检测到缓存文件：%s\n", gene_sites_rds))
  cat("直接读取 gene_sites 缓存。\n")
  
  gene_sites <- readRDS(gene_sites_rds)
  
} else {
  
  cat("\n未检测到 gene_sites 缓存，开始解析 GenBank。\n")
  
  gb_lines <- read_lines(gb_file)
  
  origin_start <- which(str_detect(gb_lines, "^ORIGIN"))
  end_idx      <- which(str_detect(gb_lines, "^//"))
  
  if (length(origin_start) == 0) stop("GenBank 文件中未找到 ORIGIN。")
  if (length(end_idx) == 0) stop("GenBank 文件中未找到 // 结束标记。")
  
  seq_lines <- gb_lines[(origin_start + 1):(end_idx[1] - 1)]
  
  genome_chr <- seq_lines |>
    str_replace_all("[0-9]", "") |>
    str_replace_all("\\s", "") |>
    paste0(collapse = "") |>
    str_to_upper()
  
  genome <- DNAString(genome_chr)
  
  cat(sprintf("基因组长度: %d bp\n", length(genome)))
  
  feat_start <- which(str_detect(gb_lines, "^FEATURES"))
  if (length(feat_start) == 0) stop("GenBank 文件中未找到 FEATURES。")
  
  feat_block <- gb_lines[(feat_start + 1):(origin_start - 1)]
  
  is_feat_head <- str_detect(feat_block, "^ {5}\\S")
  head_idx     <- which(is_feat_head)
  
  if (length(head_idx) == 0) stop("FEATURES 区块中未检测到 feature。")
  
  parse_one_feature <- function(start_i, end_i) {
    
    block <- feat_block[start_i:end_i]
    keyline <- block[1]
    ftype   <- str_trim(str_sub(keyline, 6, 21))
    
    loc_lines_raw  <- str_sub(block, 22)
    loc_lines_trim <- str_trim(loc_lines_raw)
    
    qual_idx <- which(str_detect(loc_lines_trim, "^/"))
    
    if (length(qual_idx) > 0) {
      if (qual_idx[1] == 1) {
        loc_lines <- character(0)
      } else {
        loc_lines <- loc_lines_trim[1:(qual_idx[1] - 1)]
      }
    } else {
      loc_lines <- loc_lines_trim
    }
    
    loc_lines <- loc_lines[loc_lines != ""]
    loc <- paste0(loc_lines, collapse = "")
    
    joined <- paste(block, collapse = " ")
    
    locus <- str_match(joined, '/locus_tag="([^"]+)"')[, 2]
    if (is.na(locus)) {
      locus <- str_match(joined, '/gene="([^"]+)"')[, 2]
    }
    
    gene_id <- str_extract(locus, "Rv\\d{4}[A-Za-z]?")
    
    tibble(
      type    = ftype,
      loc     = loc,
      gene_id = gene_id
    )
  }
  
  feat_tbl <- map2_dfr(
    head_idx,
    c(head_idx[-1] - 1, length(feat_block)),
    parse_one_feature
  ) |>
    filter(type == "CDS", !is.na(gene_id), !is.na(loc), loc != "")
  
  cat(sprintf("解析到 CDS 数: %d\n", nrow(feat_tbl)))
  
  extract_cds_seq <- function(loc) {
    
    if (is.na(loc) || loc == "") return(NA_character_)
    
    is_comp <- str_detect(loc, "complement")
    nums <- str_extract_all(loc, "\\d+")[[1]] |> as.integer()
    
    if (length(nums) < 2) return(NA_character_)
    if (length(nums) %% 2 != 0) {
      warning(sprintf("坐标数字数量不是偶数，跳过：%s", loc))
      return(NA_character_)
    }
    
    pairs <- matrix(nums, ncol = 2, byrow = TRUE)
    
    if (any(pairs[, 1] < 1) || any(pairs[, 2] > length(genome))) {
      warning(sprintf("坐标超出基因组范围，跳过：%s", loc))
      return(NA_character_)
    }
    
    segs <- apply(
      pairs, 1,
      function(p) as.character(subseq(genome, start = p[1], end = p[2]))
    )
    
    full <- paste0(segs, collapse = "")
    
    if (is_comp) {
      full <- as.character(reverseComplement(DNAString(full)))
    }
    
    full
  }
  
  cds_tbl <- feat_tbl |>
    mutate(
      cds_seq = map_chr(loc, extract_cds_seq),
      cds_len = nchar(cds_seq)
    ) |>
    filter(!is.na(cds_seq), cds_len >= 3) |>
    arrange(gene_id, desc(cds_len)) |>
    group_by(gene_id) |>
    slice(1) |>
    ungroup() |>
    select(gene_id, loc, cds_seq, cds_len)
  
  cat(sprintf("去重后 CDS 基因数: %d\n", nrow(cds_tbl)))
  
  bases <- c("A", "C", "G", "T")
  
  all_codons <- as.vector(
    outer(
      as.vector(outer(bases, bases, paste0)),
      bases,
      paste0
    )
  )
  
  genetic_code <- Biostrings::GENETIC_CODE
  
  make_codon_site_table <- function() {
    
    codon_tbl <- lapply(all_codons, function(codon) {
      
      aa0 <- unname(genetic_code[codon])
      syn <- 0
      nonsyn <- 0
      
      for (pos in 1:3) {
        orig <- substr(codon, pos, pos)
        
        for (b in setdiff(bases, orig)) {
          codon_mut <- codon
          substr(codon_mut, pos, pos) <- b
          aa1 <- unname(genetic_code[codon_mut])
          
          if (aa1 == aa0) {
            syn <- syn + 1
          } else {
            nonsyn <- nonsyn + 1
          }
        }
      }
      
      tibble(
        codon   = codon,
        N_sites = nonsyn / 3,
        S_sites = syn / 3
      )
    })
    
    bind_rows(codon_tbl)
  }
  
  codon_site_tbl <- make_codon_site_table()
  
  N_site_lookup <- setNames(codon_site_tbl$N_sites, codon_site_tbl$codon)
  S_site_lookup <- setNames(codon_site_tbl$S_sites, codon_site_tbl$codon)
  
  stop_codons <- c("TAA", "TAG", "TGA")
  
  count_sites_cds_fast <- function(dna_seq) {
    
    dna_seq <- toupper(dna_seq)
    L <- nchar(dna_seq)
    L3 <- L - (L %% 3)
    
    if (L3 < 3) {
      return(c(N_sites = NA_real_, S_sites = NA_real_, L_bp = L))
    }
    
    dna_seq <- substr(dna_seq, 1, L3)
    
    codons <- substring(dna_seq, seq(1, L3, 3), seq(3, L3, 3))
    
    if (length(codons) > 0 && tail(codons, 1) %in% stop_codons) {
      codons <- head(codons, -1)
    }
    
    if (length(codons) == 0) {
      return(c(N_sites = NA_real_, S_sites = NA_real_, L_bp = L3))
    }
    
    valid <- grepl("^[ACGT]{3}$", codons)
    codons <- codons[valid]
    
    if (length(codons) == 0) {
      return(c(N_sites = NA_real_, S_sites = NA_real_, L_bp = L3))
    }
    
    c(
      N_sites = sum(N_site_lookup[codons], na.rm = TRUE),
      S_sites = sum(S_site_lookup[codons], na.rm = TRUE),
      L_bp    = length(codons) * 3
    )
  }
  
  gene_sites <- cds_tbl |>
    mutate(
      sites   = map(cds_seq, count_sites_cds_fast),
      N_sites = map_dbl(sites, ~ .x[["N_sites"]]),
      S_sites = map_dbl(sites, ~ .x[["S_sites"]]),
      Li_bp   = map_dbl(sites, ~ .x[["L_bp"]])
    ) |>
    select(gene_id, N_sites, S_sites, Li_bp) |>
    filter(
      !is.na(N_sites), !is.na(S_sites), !is.na(Li_bp),
      N_sites > 0, S_sites > 0, Li_bp > 0
    )
  
  cat(sprintf("成功计算 N/S sites 的基因数: %d\n", nrow(gene_sites)))
  
  saveRDS(gene_sites, gene_sites_rds)
  cat(sprintf("gene_sites 已保存到：%s\n", gene_sites_rds))
}

# ============================================================
# 5) 基础汇总
# ============================================================

cat("\n---- gene_sites 概览 ----\n")
print(summary(gene_sites))

cat(sprintf(
  "总 N_sites = %.2f\n总 S_sites = %.2f\n总 CDS bp  = %.0f\n",
  sum(gene_sites$N_sites, na.rm = TRUE),
  sum(gene_sites$S_sites, na.rm = TRUE),
  sum(gene_sites$Li_bp,   na.rm = TRUE)
))

# ============================================================
# 6) 样本级 dN/dS
# ============================================================

mut_coding <- mut |>
  filter(!is.na(gene_id)) |>
  filter(is_nonsyn | is_syn)

cat("\n---- coding mutation 检查 ----\n")
cat(sprintf(
  "coding N = %d\ncoding S = %d\n",
  sum(mut_coding$is_nonsyn, na.rm = TRUE),
  sum(mut_coding$is_syn,    na.rm = TRUE)
))

N_sites_total_genome <- sum(gene_sites$N_sites, na.rm = TRUE)
S_sites_total_genome <- sum(gene_sites$S_sites, na.rm = TRUE)

sample_dnds <- mut_coding |>
  summarise(
    N_count = sum(is_nonsyn, na.rm = TRUE),
    S_count = sum(is_syn,    na.rm = TRUE),
    .by = sample
  ) |>
  mutate(
    N_sites_total = N_sites_total_genome,
    S_sites_total = S_sites_total_genome,
    dN = (N_count + 0.5) / N_sites_total,
    dS = (S_count + 0.5) / S_sites_total,
    dNdS = dN / dS
  ) |>
  arrange(desc(dNdS))

write_tsv(sample_dnds, "per_sample_unfixed_dNdS.tsv")

# ============================================================
# 6b) 样本级：只用发生突变的基因作为 denominator
# ============================================================

# 1. 先统计每个样本全部的 N_count 和 S_count（不做任何过滤）
all_mutation_counts <- mut_coding |>
  summarise(
    N_count = sum(is_nonsyn),
    S_count = sum(is_syn),
    .by = sample
  )

# 2. 计算只使用“发生突变的基因”的 N_sites 和 S_sites 总和
mutated_gene_sites <- mut_coding |>
  distinct(sample, gene_id) |>
  left_join(gene_sites, by = "gene_id") |>
  filter(!is.na(N_sites), !is.na(S_sites), N_sites > 0, S_sites > 0) |>
  summarise(
    N_sites_total = sum(N_sites),
    S_sites_total = sum(S_sites),
    .by = sample
  )

# 3. 合并得到最终结果（N_count 和 S_count 与 sample_dNdS 一致）
sample_dnds_mutated_genes_denominator <- all_mutation_counts |>
  left_join(mutated_gene_sites, by = "sample") |>
  mutate(
    dN = N_count / N_sites_total,
    dS = S_count / S_sites_total,
    dNdS_mut_den = dN / dS
  ) |>
  select(sample, N_count, S_count, N_sites_total, S_sites_total, dN, dS, dNdS_mut_den)

# 输出到 Excel（保留你原逻辑，并新增 sheet：sample_nonfixed_variant_counts）
writexl::write_xlsx(
  list(
    sample_dNdS = sample_dnds,
    sample_dNdS_mut_den = sample_dnds_mutated_genes_denominator,
    sample_nonfixed_variant_counts = per_sample_nonfixed_variant_counts
  ),
  path = excel_file
)

# ============================================================
# 7) 全基因组基因 burden Di（已修复：使用标准突变密度 Di = Mi / Li_bp，并添加 Poisson burden 检验）
# ============================================================

# ============================================================
# 7) 全基因组基因 burden Di（严格按照参考文献：Di = Mi / (N * Li_bp)）
# ============================================================

gene_mut_counts_for_Di <- mut_coding |>
  filter(!is.na(gene_id)) |>
  group_by(gene_id) |>
  summarise(Mi = n(), .groups = "drop")

# 计算全基因组总突变数 N
N_total_mut <- sum(gene_mut_counts_for_Di$Mi, na.rm = TRUE)

gene_burden <- gene_mut_counts_for_Di |>
  left_join(gene_sites, by = "gene_id") |>
  filter(!is.na(Li_bp), Li_bp > 0) |>
  mutate(
    Di = Mi / (N_total_mut * Li_bp),  # 严格按照参考文献公式
    Di_abs = Mi / Li_bp,             # 保留绝对密度供参考
    Mi_per_kb = Mi / (Li_bp / 1000),
    expected_rate = mean(Mi, na.rm = TRUE) / mean(Li_bp, na.rm = TRUE),
    p_burden = pmap_dbl(
      list(Mi, Li_bp),
      ~ poisson.test(..1,
                     T = ..2 / mean(gene_sites$Li_bp, na.rm = TRUE) * mean(gene_mut_counts_for_Di$Mi, na.rm = TRUE)
      )$p.value
    )
  ) |>
  arrange(desc(Di))

write_tsv(gene_burden, "mutation_burden_gene_Di.tsv")

# ============================================================
# 7b) 全基因组每个基因 dN/dS + 事件级 bootstrap CI + 双侧p值
# ============================================================

gene_mut_counts <- mut_coding |>
  summarise(
    N_count = sum(is_nonsyn, na.rm = TRUE),
    S_count = sum(is_syn,    na.rm = TRUE),
    Mi      = n(),
    .by = gene_id
  )

gene_dnds_base <- gene_sites |>
  left_join(gene_mut_counts, by = "gene_id") |>
  mutate(
    N_count = replace_na(N_count, 0),
    S_count = replace_na(S_count, 0),
    Mi      = replace_na(Mi, 0),
    
    dN      = (N_count + 0.5) / N_sites,
    dS      = (S_count + 0.5) / S_sites,
    dNdS    = dN / dS
  )

N_total_unfixed <- sum(gene_mut_counts$Mi, na.rm = TRUE)

gene_dnds_base <- gene_dnds_base |>
  mutate(
    Di = Mi / Li_bp
  )

set.seed(123)
bootstrap_reps <- 1000

events <- mut_coding |>
  filter(!is.na(gene_id)) |>
  select(gene_id, is_nonsyn, is_syn) |>
  group_by(gene_id) |>
  group_split()

gene_ids_list <- map_chr(events, ~ unique(.x$gene_id)[1])

gene_boot_rows <- vector("list", length(gene_ids_list))
names(gene_boot_rows) <- gene_ids_list

for (i in seq_along(events)) {
  g_id <- gene_ids_list[i]
  dfg  <- events[[i]]
  
  gene_row <- gene_sites |> filter(gene_id == g_id)
  if (nrow(gene_row) != 1) {
    gene_boot_rows[[i]] <- tibble(gene_id = g_id, replicate = integer(0), dNdS_boot = numeric(0))
    next
  }
  
  N_sites_g <- gene_row$N_sites[1]
  S_sites_g <- gene_row$S_sites[1]
  
  n_events <- nrow(dfg)
  if (n_events == 0 || !is.finite(N_sites_g) || !is.finite(S_sites_g) || N_sites_g <= 0 || S_sites_g <= 0) {
    gene_boot_rows[[i]] <- tibble(gene_id = g_id, replicate = integer(0), dNdS_boot = numeric(0))
    next
  }
  
  idx_list <- replicate(
    bootstrap_reps,
    sample.int(n_events, size = n_events, replace = TRUE),
    simplify = FALSE
  )
  
  dNdS_boot_vec <- vapply(idx_list, function(idx) {
    boot_dat <- dfg[idx, , drop = FALSE]
    N_count_boot <- sum(boot_dat$is_nonsyn, na.rm = TRUE)
    S_count_boot <- sum(boot_dat$is_syn,    na.rm = TRUE)
    dN_boot <- (N_count_boot + 0.5) / N_sites_g
    dS_boot <- (S_count_boot + 0.5) / S_sites_g
    dN_boot / dS_boot
  }, FUN.VALUE = numeric(1))
  
  dNdS_boot_vec <- dNdS_boot_vec[is.finite(dNdS_boot_vec) & !is.na(dNdS_boot_vec)]
  
  gene_boot_rows[[i]] <- tibble(
    gene_id = g_id,
    replicate = seq_along(dNdS_boot_vec),
    dNdS_boot = as.numeric(dNdS_boot_vec)
  )
}

gene_boot_df <- bind_rows(gene_boot_rows)

gene_boot_ci <- gene_boot_df |>
  group_by(gene_id) |>
  summarise(
    dNdS_boot_low  = ifelse(n() >= 2, as.numeric(quantile(dNdS_boot, 0.025, na.rm = TRUE)), NA_real_),
    dNdS_boot_high = ifelse(n() >= 2, as.numeric(quantile(dNdS_boot, 0.975, na.rm = TRUE)), NA_real_),
    dNdS_boot_mean = mean(dNdS_boot, na.rm = TRUE),
    .groups = "drop"
  )

gene_boot_pvals <- gene_boot_df |>
  group_by(gene_id) |>
  summarise(
    p_upper = mean(dNdS_boot >= 1, na.rm = TRUE),
    p_lower = mean(dNdS_boot <= 1, na.rm = TRUE),
    p_value = 2 * pmin(p_upper, p_lower),
    p_value = pmin(p_value, 1),
    dNdS_boot_mean = mean(dNdS_boot, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    neg_log10_p_adj = ifelse(!is.na(p_adj) & p_adj > 0, -log10(p_adj), NA_real_),
    signif_class = case_when(
      !is.na(p_adj) & p_adj < 0.05 & dNdS_boot_mean > 1 ~ "Significant (dNdS > 1)",
      !is.na(p_adj) & p_adj < 0.05 & dNdS_boot_mean < 1 ~ "Significant (dNdS < 1)",
      TRUE ~ "Not significant"
    )
  )

gene_dnds_burden <- gene_dnds_base |>
  left_join(gene_boot_ci, by = "gene_id") |>
  left_join(
    gene_boot_pvals |>
      select(gene_id, p_value, p_adj, neg_log10_p_adj, signif_class),
    by = "gene_id"
  ) |>
  mutate(
    neg_log10_p_adj = ifelse(is.na(neg_log10_p_adj), NA_real_, neg_log10_p_adj)
  ) |>
  select(
    gene_id, N_sites, S_sites, Li_bp,
    Mi, N_count, S_count,
    dN, dS, dNdS,
    dNdS_boot_low, dNdS_boot_high,
    Di,
    p_value, p_adj, neg_log10_p_adj, signif_class, dNdS_boot_mean
  ) |>
  arrange(desc(Mi))

write_tsv(gene_dnds_burden, "all_genes_dNdS_and_Di_bootCI.tsv")

# ============================================================
# 8) 候选基因集（保持你原样）
# ============================================================

candidate_genes <- tribble(
  ~module, ~gene_symbol, ~gene_id, ~function_note,
  "DosR_core", "dosR/devR", "Rv3133c", "Dormancy response regulator",
  "DosR_core", "dosS/devS", "Rv3132c", "Dormancy sensor histidine kinase",
  "DosR_core", "dosT",     "Rv2027c", "Dormancy sensor histidine kinase",
  "DosR_regulon", "hspX/acr", "Rv2031c", "Alpha-crystallin; dormancy marker",
  "DosR_regulon", "acr2",     "Rv0251c", "Stress response protein",
  "DosR_regulon", "tgs1",     "Rv3130c", "Triacylglycerol synthase; energy storage",
  "DosR_regulon", "narK2",    "Rv1737c", "Nitrate transporter; hypoxia adaptation",
  "DosR_regulon", "narX",     "Rv1736c", "Nitrate metabolism-associated protein",
  "DosR_regulon", "fdxA",     "Rv2007c", "Ferredoxin; redox metabolism",
  "DosR_regulon", "Rv2623",   "Rv2623", "DosR regulon protein",
  "DosR_regulon", "Rv2626c",  "Rv2626c", "DosR regulon protein",
  "DosR_regulon", "Rv2627c",  "Rv2627c", "DosR regulon protein",
  "DosR_regulon", "Rv2628",   "Rv2628", "DosR regulon protein",
  "DosR_regulon", "Rv1733c",  "Rv1733c", "Hypoxia-induced protein",
  "DosR_regulon", "Rv1738",   "Rv1738", "Dormancy-associated protein",
  "DosR_regulon", "Rv1813c",  "Rv1813c", "DosR regulon-associated protein",
  "DosR_regulon", "Rv1996",   "Rv1996", "DosR regulon-associated protein",
  "DosR_regulon", "Rv1997",   "Rv1997", "Hypoxia response-associated protein",
  "DosR_regulon", "Rv2003c",  "Rv2003c", "DosR regulon-associated protein",
  "DosR_regulon", "Rv2004c",  "Rv2004c", "DosR regulon-associated protein",
  "DosR_regulon", "Rv2005c",  "Rv2005c", "DosR regulon-associated protein",
  "DosR_regulon", "Rv2006",    "Rv2006", "DosR regulon-associated protein",
  "DosR_regulon", "Rv2028c",   "Rv2028c", "DosR regulon-associated protein",
  "DosR_regulon", "Rv2029c",   "Rv2029c", "DosR regulon-associated protein",
  "DosR_regulon", "Rv2030c",  "Rv2030c", "hspX-adjacent dormancy-associated protein",
  "DosR_regulon", "Rv3126c",  "Rv3126c", "DosR regulon-associated protein",
  "DosR_regulon", "Rv3127",    "Rv3127", "DosR regulon-associated protein",
  "DosR_regulon", "Rv3129",    "Rv3129", "DosR regulon-associated protein",
  "DosR_regulon", "Rv3131",    "Rv3131", "DosR regulon-associated protein",
  "Rpf_resuscitation", "rpfA", "Rv0867c", "Resuscitation-promoting factor A",
  "Rpf_resuscitation", "rpfB", "Rv1009",  "Resuscitation-promoting factor B",
  "Rpf_resuscitation", "rpfC", "Rv1884c", "Resuscitation-promoting factor C",
  "Rpf_resuscitation", "rpfD", "Rv2389c", "Resuscitation-promoting factor D",
  "Rpf_resuscitation", "rpfE", "Rv2450c", "Resuscitation-promoting factor E",
  "Cell_wall_remodeling", "ripa",  "Rv1477",  "Peptidoglycan hydrolase; interacts with RpfB",
  "Cell_wall_remodeling", "ponA1", "Rv0050",  "Penicillin-binding protein; cell wall synthesis",
  "Cell_wall_remodeling", "ponA2", "Rv3682",  "Penicillin-binding protein; cell wall synthesis",
  "Cell_wall_remodeling", "dacB1", "Rv3330",  "D-alanyl-D-alanine carboxypeptidase",
  "Cell_wall_remodeling", "dacB2", "Rv2911",   "D-alanyl-D-alanine carboxypeptidase",
  "Cell_division", "ftsZ",  "Rv2150c", "Tubulin-like cell division protein",
  "Cell_division", "ftsQ",  "Rv2151c", "Cell division protein FtsQ",
  "Cell_division", "ftsW",  "Rv2154c", "Cell division protein FtsW",
  "Cell_division", "ftsI",  "Rv2163c", "Cell division penicillin-binding protein",
  "Cell_division", "sepF",  "Rv2147c", "Cell division-associated protein SepF",
  "Peptidoglycan_synthesis", "murA", "Rv1315", "UDP-N-acetylglucosamine enolpyruvyl transferase",
  "Peptidoglycan_synthesis", "murB", "Rv0482", "UDP-N-acetylenolpyruvoylglucosamine reductase",
  "Peptidoglycan_synthesis", "murC", "Rv2152c", "UDP-N-acetylmuramate-L-alanine ligase",
  "Peptidoglycan_synthesis", "murD", "Rv2155c", "UDP-N-acetylmuramoylalanine-D-glutamate ligase",
  "Peptidoglycan_synthesis", "murE", "Rv2158c", "UDP-N-acetylmuramoylalanyl-D-glutamate-2,6-diaminopimelate ligase",
  "Peptidoglycan_synthesis", "murF", "Rv2157c", "UDP-N-acetylmuramoyl-tripeptide-D-alanyl-D-alanine ligase",
  "Peptidoglycan_synthesis", "ddlA", "Rv2981c", "D-alanine-D-alanine ligase",
  "Peptidoglycan_synthesis", "alr",  "Rv3423c", "Alanine racemase",
  "Mycolic_acid_envelope", "kasA", "Rv2245", "Beta-ketoacyl-ACP synthase",
  "Mycolic_acid_envelope", "kasB", "Rv2246", "Beta-ketoacyl-ACP synthase",
  "Mycolic_acid_envelope", "inhA", "Rv1484", "Enoyl-ACP reductase",
  "Mycolic_acid_envelope", "fabD", "Rv2243", "Malonyl CoA-acyl carrier protein transacylase",
  "Mycolic_acid_envelope", "acpM", "Rv2244", "Acyl carrier protein",
  "Mycolic_acid_envelope", "mmaA1", "Rv0645c", "Mycolic acid methyltransferase",
  "Mycolic_acid_envelope", "mmaA2", "Rv0644c", "Mycolic acid methyltransferase",
  "Mycolic_acid_envelope", "mmaA3", "Rv0643c", "Mycolic acid methyltransferase",
  "Mycolic_acid_envelope", "mmaA4", "Rv0642c", "Mycolic acid methyltransferase"
) |>
  distinct(module, gene_id, .keep_all = TRUE)

write_tsv(candidate_genes, "candidate_dormancy_resuscitation_genes.tsv")

# ============================================================
# 9) 候选基因检查
# ============================================================

candidate_gene_check <- candidate_genes |>
  left_join(gene_sites, by = "gene_id") |>
  mutate(in_gene_sites = !is.na(Li_bp)) |>
  arrange(module, gene_id)

write_tsv(candidate_gene_check, "candidate_gene_check_in_H37Rv_gene_sites.tsv")

cat("\n---- 候选基因在 gene_sites 中的匹配情况 ----\n")
print(
  candidate_gene_check |>
    group_by(module, in_gene_sites) |>
    summarise(n = n(), .groups = "drop") |>
    arrange(module, desc(in_gene_sites))
)

cat("\n---- 未匹配到 gene_sites 的候选基因 ----\n")
print(candidate_gene_check |>
        filter(!in_gene_sites) |>
        select(module, gene_symbol, gene_id, function_note))

# ============================================================
# 10) 候选模块 mutation burden（已同步修复 Di 计算）
# ============================================================

candidate_gene_burden <- candidate_genes |>
  left_join(gene_burden, by = "gene_id") |>
  left_join(
    mut_coding |>
      summarise(
        N_count = sum(is_nonsyn, na.rm = TRUE),
        S_count = sum(is_syn,    na.rm = TRUE),
        .by = gene_id
      ),
    by = "gene_id"
  ) |>
  mutate(
    Mi = replace_na(Mi, 0),
    N_count = replace_na(N_count, 0),
    S_count = replace_na(S_count, 0),
    
    Di = dplyr::if_else(
      !is.na(Li_bp) & Li_bp > 0,
      Mi / Li_bp,
      NA_real_
    ),
    Mi_per_kb = dplyr::if_else(!is.na(Li_bp) & Li_bp > 0, Mi / (Li_bp / 1000), NA_real_),
    N_per_kb = dplyr::if_else(!is.na(Li_bp) & Li_bp > 0, N_count / (Li_bp / 1000), NA_real_),
    S_per_kb = dplyr::if_else(!is.na(Li_bp) & Li_bp > 0, S_count / (Li_bp / 1000), NA_real_)
  ) |>
  arrange(module, desc(Mi_per_kb))

write_tsv(candidate_gene_burden, "candidate_gene_mutation_burden.tsv")

# ============================================================
# 11) 模块级 mutation burden（已同步修复 Di）
# ============================================================

candidate_gene_ids <- unique(candidate_genes$gene_id)

module_sites <- candidate_genes |>
  left_join(gene_sites, by = "gene_id") |>
  filter(!is.na(Li_bp), Li_bp > 0) |>
  summarise(
    module_gene_number = n_distinct(gene_id),
    module_Li_bp       = sum(Li_bp, na.rm = TRUE),
    module_N_sites     = sum(N_sites, na.rm = TRUE),
    module_S_sites     = sum(S_sites, na.rm = TRUE),
    .by = module
  )

module_mut_counts <- mut |>
  filter(!is.na(gene_id)) |>
  inner_join(candidate_genes |> select(module, gene_id), by = "gene_id") |>
  summarise(Mi = n(), .by = module)

module_NS_counts <- mut_coding |>
  inner_join(candidate_genes |> select(module, gene_id), by = "gene_id") |>
  summarise(
    N_count = sum(is_nonsyn, na.rm = TRUE),
    S_count = sum(is_syn,    na.rm = TRUE),
    .by = module
  )

module_burden <- module_sites |>
  left_join(module_mut_counts, by = "module") |>
  left_join(module_NS_counts, by = "module") |>
  mutate(
    Mi = replace_na(Mi, 0),
    N_count = replace_na(N_count, 0),
    S_count = replace_na(S_count, 0),
    
    module_Mi_per_kb = dplyr::if_else(
      !is.na(module_Li_bp) & module_Li_bp > 0,
      Mi / (module_Li_bp / 1000),
      NA_real_
    ),
    
    module_Di = dplyr::if_else(
      !is.na(module_Li_bp) & module_Li_bp > 0,
      Mi / module_Li_bp,
      NA_real_
    ),
    
    dN = (N_count + 0.5) / module_N_sites,
    dS = (S_count + 0.5) / module_S_sites,
    dNdS = dN / dS
  ) |>
  arrange(desc(module_Mi_per_kb))

write_tsv(module_burden, "candidate_module_mutation_burden.tsv")

# ============================================================
# 12-13) Fisher 检验（保持你原逻辑）
# ============================================================

background_coding <- mut_coding |>
  filter(!gene_id %in% candidate_gene_ids)

background_N <- sum(background_coding$is_nonsyn, na.rm = TRUE)
background_S <- sum(background_coding$is_syn,    na.rm = TRUE)

fisher_one_module <- function(module_name) {
  
  x <- mut_coding |>
    inner_join(
      candidate_genes |> filter(module == module_name) |> select(module, gene_id),
      by = "gene_id"
    )
  
  module_N <- sum(x$is_nonsyn, na.rm = TRUE)
  module_S <- sum(x$is_syn,    na.rm = TRUE)
  
  tab <- matrix(
    c(module_N, module_S, background_N, background_S),
    nrow = 2,
    byrow = TRUE
  )
  
  colnames(tab) <- c("Nonsynonymous", "Synonymous")
  rownames(tab) <- c(module_name, "Background")
  
  if (sum(tab[1, ]) == 0 || sum(tab[2, ]) == 0) {
    return(tibble(
      module = module_name,
      module_N = module_N,
      module_S = module_S,
      background_N = background_N,
      background_S = background_S,
      odds_ratio = NA_real_,
      p_value = NA_real_
    ))
  }
  
  ft <- fisher.test(tab)
  
  tibble(
    module = module_name,
    module_N = module_N,
    module_S = module_S,
    background_N = background_N,
    background_S = background_S,
    odds_ratio = as.numeric(ft$estimate),
    p_value = ft$p.value
  )
}

module_fisher <- map_dfr(unique(candidate_genes$module), fisher_one_module) |>
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) |>
  arrange(p_value)

write_tsv(module_fisher, "candidate_module_nonsyn_enrichment_fisher.tsv")

candidate_coding <- mut_coding |>
  filter(gene_id %in% candidate_gene_ids)

candidate_all_N <- sum(candidate_coding$is_nonsyn, na.rm = TRUE)
candidate_all_S <- sum(candidate_coding$is_syn,    na.rm = TRUE)

candidate_vs_background_tab <- matrix(
  c(candidate_all_N, candidate_all_S, background_N, background_S),
  nrow = 2,
  byrow = TRUE
)

colnames(candidate_vs_background_tab) <- c("Nonsynonymous", "Synonymous")
rownames(candidate_vs_background_tab) <- c("Candidate_genes", "Background")

candidate_vs_background_fisher <- fisher.test(candidate_vs_background_tab)

candidate_vs_background_result <- tibble(
  comparison = "Candidate dormancy-resuscitation genes vs background",
  candidate_N = candidate_all_N,
  candidate_S = candidate_all_S,
  background_N = background_N,
  background_S = background_S,
  odds_ratio = as.numeric(candidate_vs_background_fisher$estimate),
  p_value = candidate_vs_background_fisher$p.value
)

write_tsv(candidate_vs_background_result, "candidate_genes_vs_background_nonsyn_enrichment.tsv")

# ============================================================
# 14) 样本级：候选基因集 mutation burden 和 dN/dS（保持你原逻辑）
# ============================================================

sample_candidate_summary <- mut_coding |>
  mutate(in_candidate = gene_id %in% candidate_gene_ids) |>
  summarise(
    candidate_Mi = sum(in_candidate, na.rm = TRUE),
    candidate_N  = sum(in_candidate & is_nonsyn, na.rm = TRUE),
    candidate_S  = sum(in_candidate & is_syn,    na.rm = TRUE),
    background_Mi = sum(!in_candidate, na.rm = TRUE),
    background_N  = sum(!in_candidate & is_nonsyn, na.rm = TRUE),
    background_S  = sum(!in_candidate & is_syn,    na.rm = TRUE),
    .by = sample
  )

candidate_sites_total <- candidate_genes |>
  left_join(gene_sites, by = "gene_id") |>
  filter(!is.na(Li_bp), Li_bp > 0) |>
  summarise(
    candidate_Li_bp   = sum(Li_bp, na.rm = TRUE),
    candidate_N_sites = sum(N_sites, na.rm = TRUE),
    candidate_S_sites = sum(S_sites, na.rm = TRUE)
  )

background_sites_total <- gene_sites |>
  filter(!gene_id %in% candidate_gene_ids) |>
  summarise(
    background_Li_bp   = sum(Li_bp, na.rm = TRUE),
    background_N_sites = sum(N_sites, na.rm = TRUE),
    background_S_sites = sum(S_sites, na.rm = TRUE)
  )

sample_candidate_summary <- sample_candidate_summary |>
  mutate(
    candidate_Li_bp   = candidate_sites_total$candidate_Li_bp,
    candidate_N_sites = candidate_sites_total$candidate_N_sites,
    candidate_S_sites = candidate_sites_total$candidate_S_sites,
    background_Li_bp   = background_sites_total$background_Li_bp,
    background_N_sites = background_sites_total$background_N_sites,
    background_S_sites = background_sites_total$background_S_sites,
    candidate_Mi_per_kb = candidate_Mi / (candidate_Li_bp / 1000),
    background_Mi_per_kb = background_Mi / (background_Li_bp / 1000),
    candidate_dN = (candidate_N + 0.5) / candidate_N_sites,
    candidate_dS = (candidate_S + 0.5) / candidate_S_sites,
    candidate_dNdS = candidate_dN / candidate_dS,
    background_dN = (background_N + 0.5) / background_N_sites,
    background_dS = (background_S + 0.5) / background_S_sites,
    background_dNdS = background_dN / background_dS
  ) |>
  arrange(desc(candidate_Mi_per_kb))

write_tsv(sample_candidate_summary, "per_sample_candidate_genes_mutation_burden_dNdS.tsv")

# ============================================================
# 15) 写入 Excel 多个 sheet
# ============================================================

if (!exists("excel_file")) {
  excel_file <- "unfixed_mutation_analysis_results.xlsx"
}

wb <- openxlsx::createWorkbook()

addWorksheet(wb, "sample_dNdS")
openxlsx::writeData(wb, "sample_dNdS", sample_dnds)

addWorksheet(wb, "sample_dNdS_mut_den")
openxlsx::writeData(wb, "sample_dNdS_mut_den", sample_dnds_mutated_genes_denominator)

# 新增：写入每个样本非固定性突变数量
addWorksheet(wb, "sample_nonfixed_variant_counts")
openxlsx::writeData(wb, "sample_nonfixed_variant_counts", per_sample_nonfixed_variant_counts)

addWorksheet(wb, "gene_burden")
openxlsx::writeData(wb, "gene_burden", gene_burden)

addWorksheet(wb, "gene_dNdS_Di_bootCI")
openxlsx::writeData(wb, "gene_dNdS_Di_bootCI", gene_dnds_burden)

addWorksheet(wb, "candidate_genes")
openxlsx::writeData(wb, "candidate_genes", candidate_genes)

addWorksheet(wb, "candidate_gene_check")
openxlsx::writeData(wb, "candidate_gene_check", candidate_gene_check)

addWorksheet(wb, "candidate_gene_burden")
openxlsx::writeData(wb, "candidate_gene_burden", candidate_gene_burden)

addWorksheet(wb, "module_burden")
openxlsx::writeData(wb, "module_burden", module_burden)

addWorksheet(wb, "module_fisher")
openxlsx::writeData(wb, "module_fisher", module_fisher)

addWorksheet(wb, "candidate_vs_background")
openxlsx::writeData(wb, "candidate_vs_background", candidate_vs_background_result)

addWorksheet(wb, "sample_candidate_summary")
openxlsx::writeData(wb, "sample_candidate_summary", sample_candidate_summary)

addWorksheet(wb, "gene_sites")
openxlsx::writeData(wb, "gene_sites", gene_sites)

openxlsx::saveWorkbook(wb, excel_file, overwrite = TRUE)

# ============================================================
# 16) 也输出 TSV
# ============================================================

write_tsv(gene_sites, "H37Rv_gene_sites.tsv")

# ============================================================
# 17) 最终打印
# ============================================================

cat("\n完成。\n")

cat("\n---- 样本级 dN/dS，前 10 行 ----\n")
print(head(sample_dnds, 10))

cat("\n---- 样本级 dN/dS，突变基因 denominator，前 10 行 ----\n")
print(head(sample_dnds_mutated_genes_denominator, 10))

cat("\n---- 每个样本 非固定性突变数量，前 10 行 ----\n")
print(head(per_sample_nonfixed_variant_counts, 10))

cat("\n---- 全基因组 gene burden，前 10 行 ----\n")
print(head(gene_burden, 10))

cat("\n---- 全基因组每个基因 dN/dS + Di + bootstrap CI，前 10 行 ----\n")
print(head(gene_dnds_burden, 10))

cat("\n---- 候选基因 burden，前 20 行 ----\n")
print(head(candidate_gene_burden, 20))

cat("\n---- 候选模块 burden ----\n")
print(module_burden)

cat("\n---- 候选模块 nonsynonymous 富集检验 ----\n")
print(module_fisher)

cat("\n---- 候选基因集 vs 背景基因 Fisher 检验 ----\n")
print(candidate_vs_background_result)

cat(sprintf("\nExcel 已输出：%s\n", excel_file))


# ============================================================
# 18) 全基因组 gene-level nonsynonymous enrichment volcano（严格检验 dN/dS = 1，保留全部基因）
# ============================================================

# 0. 读取 Tuberculist 功能注释文件
tuberculist_anno <- read_tsv(
  "/home/data/t040512/UnfixedMutation/Tuberculist_new_20150307.txt",
  col_names = c("gene_id", "gene_symbol_tb", "strand", "start", "end", "function_note", "category"),
  col_types = cols(
    gene_id = col_character(),
    gene_symbol_tb = col_character(),
    strand = col_character(),
    start = col_integer(),
    end = col_integer(),
    function_note = col_character(),
    category = col_character()
  )
) |>
  # 清理 category 列的方括号，使其更整洁 (例如 "[3-Information pathways]" -> "3-Information pathways")
  mutate(category = str_replace(category, "^\\[(.*)\\]$", "\\1")) |>
  select(gene_id, gene_symbol_tb, function_note, category)

# 1. 获取全基因组所有基因 ID (保留全部基因)
all_gene_ids <- gene_sites$gene_id |> unique()

# 2. 统计每个基因的突变计数并与基因位点信息合并
gene_ns_all <- mut_coding |>
  filter(!is.na(gene_id)) |>
  summarise(
    gene_N = sum(is_nonsyn, na.rm = TRUE),
    gene_S = sum(is_syn,    na.rm = TRUE),
    .by = gene_id
  ) |>
  right_join(tibble(gene_id = all_gene_ids), by = "gene_id") |> # right_join 确保保留所有基因
  mutate(
    gene_N = replace_na(gene_N, 0),
    gene_S = replace_na(gene_S, 0)
  ) |>
  left_join(gene_sites, by = "gene_id")

# 3. 自定义函数：在对数尺度下计算二项精确检验的双侧 p 值
calc_binom_p <- function(x, n, p) {
  if (n == 0) return(NA)       # 没有突变的基因无法检验
  if (is.na(p) || p == 0 || p == 1) return(NA) # 位点信息缺失或极端情况无法检验
  
  # 在对数尺度下计算 P(X <= x) 和 P(X >= x)
  log_p_lower <- pbinom(x, n, p, lower.tail = TRUE, log.p = TRUE)
  log_p_upper <- pbinom(x - 1, n, p, lower.tail = FALSE, log.p = TRUE)
  
  # 双侧 p 值近似：2 * min(p_lower, p_upper)
  log_p_value <- log(2) + min(log_p_lower, log_p_upper)
  
  # p 值最大不超过 1 (log(1) = 0)
  if (log_p_value > 0) log_p_value <- 0
  
  # 将 log p 转换回正常 p 值
  p_val <- exp(log_p_value)
  
  # 如果依然下溢为 0，则将其设置为 R 能表示的最小双精度浮点数
  if (p_val == 0) {
    p_val <- .Machine$double.xmin  # 约等于 2.2e-308
  }
  
  return(p_val)
}

# 4. 逐基因检验 dN/dS 是否偏离 1
cat("\n[进度] 开始逐基因二项精确检验 (检验 dN/dS = 1)...\n")

gene_fisher_results <- gene_ns_all |>
  rowwise() |>
  mutate(
    # 判断该基因是否具备计算条件 (有突变且有位点注释)
    can_calc = (gene_N + gene_S > 0) & !is.na(N_sites) & !is.na(S_sites) & (N_sites + S_sites > 0),
    
    # 该基因自身的期望非同义突变概率 (基于位点归一化)
    p_expected_gene = ifelse(can_calc, (N_sites + 0.5) / (N_sites + S_sites + 1), NA),
    
    # 使用自定义函数计算 p 值
    p_value = ifelse(can_calc, calc_binom_p(gene_N, gene_N + gene_S, p_expected_gene), NA),
    
    # 基于位点归一化的 dN, dS 计算
    dN = ifelse(can_calc, (gene_N + 0.5) / (N_sites + 0.5), NA),
    dS = ifelse(can_calc, (gene_S + 0.5) / (S_sites + 0.5), NA),
    
    # 直接计算 dN/dS 比值
    dNdS_ratio = ifelse(can_calc, dN / dS, NA),
    
    # 计算 log2 转换值
    log2FC_dnds = ifelse(can_calc, log2(dNdS_ratio), NA)
  ) |>
  ungroup() |>
  # 合并 Tuberculist 注释
  left_join(tuberculist_anno, by = "gene_id") |>
  # 合并自定义模块 (如果 candidate_genes 中有 module 信息)
  left_join(
    candidate_genes |> select(gene_id, module),
    by = "gene_id"
  ) |>
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    neg_log10_p = -log10(p_adj_BH),
    
    # 判定标准：结合 FDR < 0.05 和 dN/dS 与 1 的比较
    signif = case_when(
      is.na(p_adj_BH) ~ "No mutations/No sites", # 缺乏计算条件的基因
      p_adj_BH < 0.05 & dNdS_ratio < 1 ~ "Purifying selection",
      p_adj_BH < 0.05 & dNdS_ratio > 1 ~ "Positive selection",
      TRUE ~ "Neutral evolution" # FDR >= 0.05
    ),
    plot_group = ifelse(is.na(module), "Other genes", module)
  ) |>
  # 整理输出列顺序 (将注释信息放在前面，方便查看)
  select(
    gene_id, gene_symbol_tb, function_note, category, module,
    gene_N, gene_S, N_sites, S_sites, 
    any_of("Li_bp"), 
    p_expected_gene, p_value, dN, dS, dNdS_ratio, log2FC_dnds, p_adj_BH, neg_log10_p, signif,
    plot_group
  ) |>
  arrange(p_adj_BH)

# 5. 输出 Excel 与 TSV
gene_fisher_xlsx <- "gene_fisher_results.xlsx"

wb2 <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb2, "gene_fisher_results")
openxlsx::writeData(wb2, "gene_fisher_results", gene_fisher_results)
openxlsx::addWorksheet(wb2, "candidate_genes")
openxlsx::writeData(wb2, "candidate_genes", candidate_genes)
openxlsx::addWorksheet(wb2, "gene_sites")
openxlsx::writeData(wb2, "gene_sites", gene_sites)
openxlsx::saveWorkbook(wb2, gene_fisher_xlsx, overwrite = TRUE)

write_tsv(gene_fisher_results, "gene_fisher_results_normalized.tsv")

cat(sprintf("\n已输出：%s (共 %d 个基因)\n", gene_fisher_xlsx, nrow(gene_fisher_results)))

# 6. 绘制火山图 (自动过滤掉无法计算的 NA 基因)
plot_df <- gene_fisher_results |>
  filter(!is.na(log2FC_dnds), !is.na(neg_log10_p))

module_colors <- c(
  "DosR_core" = "#E41A1C",
  "DosR_regulon" = "#377EB8",
  "Rpf_resuscitation" = "#4DAF4A",
  "Cell_wall_remodeling" = "#984EA3",
  "Cell_division" = "#FF7F00",
  "Peptidoglycan_synthesis" = "#A65628",
  "Mycolic_acid_envelope" = "#F781BF",
  "Other genes" = "grey70"
)

missing_groups <- setdiff(unique(plot_df$plot_group), names(module_colors))
if (length(missing_groups) > 0) {
  module_colors <- c(module_colors, setNames(rainbow(length(missing_groups)), missing_groups))
}

p_gene_volcano <- ggplot(plot_df, aes(x = log2FC_dnds, y = neg_log10_p, color = plot_group)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = module_colors) +
  theme_bw(base_size = 12) +
  labs(
    x = expression(log[2]~(dN/dS)~"(Normalized)"),
    y = "-log10(FDR)",
    color = "Module / gene set",
    title = "Genome-wide gene-level dN/dS volcano plot"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(p_gene_volcano)

ggsave(
  filename = "gene_nonsyn_enrichment_volcano_normalized.png",
  plot = p_gene_volcano,
  width = 10,
  height = 7,
  dpi = 300
)

cat("\n---- gene_fisher_results 前 10 行 ----\n")
print(head(gene_fisher_results, 10))

cat("\n[完成] 第 18 节全基因组基因水平 dN/dS 检验已执行完毕。\n")


# ============================================================
# 19) dNdS vs 1 双侧检验 + 火山图（基因级 + 样本级）
# ============================================================

# 19.1 基因级：双侧 p_value（基于 dNdS_boot 相对 1）+ BH
if (!exists("gene_boot_df") || !exists("gene_dnds_burden")) {
  stop("gene_boot_df 或 gene_dnds_burden 未找到，请确认 7b) 已成功运行。")
}

gene_boot_df2 <- gene_boot_df |>
  filter(!is.na(dNdS_boot), is.finite(dNdS_boot))

gene_boot_pvals <- gene_boot_df2 |>
  group_by(gene_id) |>
  summarise(
    p_upper = mean(dNdS_boot >= 1, na.rm = TRUE),
    p_lower = mean(dNdS_boot <= 1, na.rm = TRUE),
    p_value = 2 * pmin(p_upper, p_lower),
    p_value = pmin(p_value, 1),
    dNdS_boot_mean = mean(dNdS_boot, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    neg_log10_p_adj = ifelse(!is.na(p_adj) & p_adj > 0, -log10(p_adj), NA_real_),
    signif_class = case_when(
      !is.na(p_adj) & p_adj < 0.05 & dNdS_boot_mean > 1 ~ "Significant (dNdS > 1)",
      !is.na(p_adj) & p_adj < 0.05 & dNdS_boot_mean < 1 ~ "Significant (dNdS < 1)",
      TRUE ~ "Not significant"
    )
  )

gene_plot_df <- gene_dnds_burden |>
  select(gene_id, dNdS, dNdS_boot_low, dNdS_boot_high) |>
  left_join(
    gene_boot_pvals |>
      select(gene_id, p_value, p_adj, neg_log10_p_adj, signif_class, dNdS_boot_mean),
    by = "gene_id"
  ) |>
  arrange(gene_id)

write_tsv(gene_plot_df, "gene_dNdS_vs1_pvalues_double_sided.tsv")

colors_gene <- c(
  "Significant (dNdS > 1)" = "#E41A1C",
  "Significant (dNdS < 1)" = "#377EB8",
  "Not significant" = "grey60"
)

p_gene_volcano <- ggplot(gene_plot_df, aes(x = dNdS, y = neg_log10_p_adj, color = signif_class)) +
  geom_point(alpha = 0.75, size = 1.8, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "solid", color = "grey50", linewidth = 0.4) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = colors_gene) +
  labs(
    x = "dNdS (point estimate; cutoff = 1)",
    y = "-log10(p_adj, two-sided; BH corrected)",
    color = "Gene test",
    title = "Gene-level: dNdS vs 1 (double-sided; x=dNdS; y=-log10(p_adj))"
  ) +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

print(p_gene_volcano)

ggsave("gene_volcano_dNdS_vs1_padj_double_sided.png",
       p_gene_volcano, width = 10, height = 7, dpi = 300)

# 19.2 样本级：dNdS vs 1 双侧检验 + 火山图（分母=样本 active genes 的 N_sites/S_sites）
if (!exists("gene_sites")) stop("gene_sites 未找到。")
if (!exists("mut_coding")) stop("mut_coding 未找到。")
if (!exists("sample_dnds") || is.null(sample_dnds$sample)) stop("sample_dnds 未找到或缺列 sample。")

# ============================================================
# （新增）每个样本：非固定性突变数量 + N/S + intergenic 非固定性数量
# ============================================================

# intergenic 的判定尽量稳健：优先用 is_intergenic；否则用 mutation_type 文本；否则置 0
mut_coding2 <- mut_coding

if (!("is_intergenic" %in% names(mut_coding2))) {
  if ("mutation_type" %in% names(mut_coding2)) {
    mut_coding2 <- mut_coding2 |>
      mutate(
        is_intergenic = case_when(
          str_detect(str_to_lower(as.character(mutation_type)), "intergenic") ~ TRUE,
          TRUE ~ FALSE
        )
      )
  } else {
    mut_coding2$is_intergenic <- FALSE
  }
}

# Nonsyn / Syn 已在后续作为 is_nonsyn / is_syn 使用；若不存在则报错
req_cols <- c("sample", "is_nonsyn", "is_syn", "is_intergenic")
missing_cols <- setdiff(req_cols, names(mut_coding2))
if (length(missing_cols) > 0) {
  stop(paste0("mut_coding 缺少必要列：", paste(missing_cols, collapse = ", ")))
}

per_sample_nonfixed_counts <- mut_coding2 |>
  group_by(sample) |>
  summarise(
    nonfixed_total = sum(is_nonsyn | is_syn, na.rm = TRUE),   # coding nonsyn+syn 事件总数
    nonsyn_nonfixed = sum(is_nonsyn, na.rm = TRUE),
    syn_nonfixed    = sum(is_syn,    na.rm = TRUE),
    intergenic_nonfixed = sum(is_intergenic, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(sample)

write_tsv(per_sample_nonfixed_counts, "per_sample_nonfixed_variant_counts.tsv")

cat("\n---- 已输出：per_sample_nonfixed_variant_counts.tsv ----\n")
print(head(per_sample_nonfixed_counts, 10))

# ============================================================
# 后续：dNdS 计算（保持你原逻辑）
# ============================================================

set.seed(123)
bootstrap_reps <- 1000

# 样本点估计：active genes 的 N_sites/S_sites
sample_point <- mut_coding |>
  filter(!is.na(gene_id), !is.na(sample)) |>
  group_by(sample) |>
  summarise(
    active_genes = list(unique(gene_id[!is.na(gene_id)])),
    N_count = sum(is_nonsyn, na.rm = TRUE),
    S_count = sum(is_syn,    na.rm = TRUE),
    .groups = "drop"
  ) |>
  rowwise() |>
  mutate(
    N_sites = sum(gene_sites$N_sites[match(unlist(active_genes), gene_sites$gene_id)], na.rm = TRUE),
    S_sites = sum(gene_sites$S_sites[match(unlist(active_genes), gene_sites$gene_id)], na.rm = TRUE),
    dN = (N_count + 0.5) / N_sites,
    dS = (S_count + 0.5) / S_sites,
    dNdS = dN / dS
  ) |>
  ungroup() |>
  select(sample, dNdS, N_count, S_count, N_sites, S_sites)

sample_boot <- map_dfr(sort(unique(mut_coding$sample)), function(smp) {
  
  df_one <- mut_coding |>
    filter(sample == smp) |>
    select(gene_id, is_nonsyn, is_syn)
  
  active_genes <- unique(df_one$gene_id[!is.na(df_one$gene_id)])
  
  N_sites_s <- sum(gene_sites$N_sites[match(active_genes, gene_sites$gene_id)], na.rm = TRUE)
  S_sites_s <- sum(gene_sites$S_sites[match(active_genes, gene_sites$gene_id)], na.rm = TRUE)
  
  n_events <- nrow(df_one)
  
  if (n_events == 0 || !is.finite(N_sites_s) || !is.finite(S_sites_s) || N_sites_s <= 0 || S_sites_s <= 0) {
    return(tibble(
      sample = smp,
      p_value = NA_real_,
      p_adj = NA_real_,
      neg_log10_p_adj = NA_real_,
      signif_class = NA_character_,
      dNdS_boot_mean = NA_real_
    ))
  }
  
  idx_list <- replicate(
    bootstrap_reps,
    sample.int(n_events, size = n_events, replace = TRUE),
    simplify = FALSE
  )
  
  dNdS_boot_vec <- vapply(idx_list, function(idx) {
    boot_dat <- df_one[idx, , drop = FALSE]
    N_count_boot <- sum(boot_dat$is_nonsyn, na.rm = TRUE)
    S_count_boot <- sum(boot_dat$is_syn,    na.rm = TRUE)
    
    dN_boot <- (N_count_boot + 0.5) / N_sites_s
    dS_boot <- (S_count_boot + 0.5) / S_sites_s
    dN_boot / dS_boot
  }, FUN.VALUE = numeric(1))
  
  dNdS_boot_vec <- dNdS_boot_vec[is.finite(dNdS_boot_vec)]
  
  if (length(dNdS_boot_vec) == 0) {
    return(tibble(
      sample = smp,
      p_value = NA_real_,
      p_adj = NA_real_,
      neg_log10_p_adj = NA_real_,
      signif_class = NA_character_,
      dNdS_boot_mean = NA_real_
    ))
  }
  
  p_upper <- mean(dNdS_boot_vec >= 1, na.rm = TRUE)
  p_lower <- mean(dNdS_boot_vec <= 1, na.rm = TRUE)
  
  p_value <- 2 * pmin(p_upper, p_lower)
  p_value <- min(p_value, 1)
  
  dNdS_boot_mean <- mean(dNdS_boot_vec, na.rm = TRUE)
  
  tibble(
    sample = smp,
    p_value = p_value,
    dNdS_boot_mean = dNdS_boot_mean
  )
})

sample_boot <- sample_boot |>
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    p_adj = ifelse(is.na(p_adj), NA_real_, pmax(p_adj, 1e-300)),
    neg_log10_p_adj = -log10(p_adj),
    dNdS_point = sample_point$dNdS[match(sample, sample_point$sample)],
    signif_class = case_when(
      !is.na(p_adj) & p_adj < 0.05 & dNdS_point > 1 ~ "Significant (dNdS > 1)",
      !is.na(p_adj) & p_adj < 0.05 & dNdS_point < 1 ~ "Significant (dNdS < 1)",
      TRUE ~ "Not significant"
    )
  ) |>
  select(sample, p_value, p_adj, neg_log10_p_adj, signif_class, dNdS_boot_mean)

sample_plot_df <- sample_point |>
  left_join(sample_boot, by = "sample") |>
  arrange(sample)

write_tsv(sample_plot_df, "sample_dNdS_vs1_pvalues_double_sided.tsv")

colors_sample <- c(
  "Significant (dNdS > 1)" = "#E41A1C",
  "Significant (dNdS < 1)" = "#377EB8",
  "Not significant" = "grey60"
)

p_sample_volcano <- ggplot(sample_plot_df, aes(x = dNdS, y = neg_log10_p_adj, color = signif_class)) +
  geom_point(alpha = 0.85, size = 2.5, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "solid", color = "grey50", linewidth = 0.4) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = colors_sample) +
  labs(
    x = "dNdS (point estimate; cutoff = 1)",
    y = "-log10(p_adj, two-sided; BH corrected)",
    color = "Sample test",
    title = "Sample-level: dNdS vs 1 (double-sided; x=dNdS; y=-log10(p_adj))"
  ) +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

print(p_sample_volcano)

ggsave(
  "sample_volcano_dNdS_vs1_padj_double_sided.png",
  p_sample_volcano, width = 10, height = 7, dpi = 300
)

cat("\n---- 输出文件 ----\n")
cat("- gene_dNdS_vs1_pvalues_double_sided.tsv\n")
cat("- gene_volcano_dNdS_vs1_padj_double_sided.png\n")
cat("- sample_dNdS_vs1_pvalues_double_sided.tsv\n")
cat("- sample_volcano_dNdS_vs1_padj_double_sided.png\n")
cat("- ", excel_file, "\n", sep = "")

cat("\n脚本运行结束。 \n")












# ============================================================
# 0) 加载包
# ============================================================
library(tidyr)
library(Biostrings)
library(readxl)
library(dplyr)
library(stringr)
library(janitor)
library(readr)
library(purrr)
library(ggplot2)

# 你的脚本前面应当已经在环境中准备好：
# gene_dnds_burden、mut_coding、gene_sites
stopifnot(exists("gene_dnds_burden"),
          exists("mut_coding"),
          exists("gene_sites"))

# ============================================================
# 追加：根据 all_genes_dNdS_and_Di_bootCI.tsv 的 signif_class 三组
# 做样本级 sample_dNdS（coding N/S），并输出 tsv + 火山图
# ============================================================

# 取三组 signif_class（去NA并保留唯一值）
signif_levels <- gene_dnds_burden %>%
  dplyr::distinct(signif_class) %>%
  dplyr::pull(signif_class) %>%
  unique()
signif_levels <- signif_levels[!is.na(signif_levels)]

cat("\n---- signif_class 分组类别 ----\n")
print(signif_levels)

# 函数：对某一组 genes 计算 sample_dNdS vs 1（双侧，BH校正）并输出
run_sample_dnds_vs1_for_gene_set <- function(genes_vec, signif_label,
                                             bootstrap_reps = 1000,
                                             seed = 123) {
  
  # 仅 coding 事件（mut_coding 中应已有 is_nonsyn / is_syn）
  mut_sub <- mut_coding %>%
    dplyr::filter(!is.na(gene_id), gene_id %in% genes_vec)
  
  if (nrow(mut_sub) == 0) {
    warning(sprintf("signif_class=%s：该组 genes 在 mut_coding 中没有任何事件，跳过。", signif_label))
    return(NULL)
  }
  
  # 样本列表
  samples_all <- sort(unique(mut_sub$sample))
  
  # ---------- 1) 点估计：样本内 dN/dS ----------
  sample_point_sub <- mut_sub %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      active_genes = list(unique(gene_id[!is.na(gene_id)])),
      N_count = sum(is_nonsyn, na.rm = TRUE),
      S_count = sum(is_syn,    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    dplyr::mutate(
      N_sites = sum(
        gene_sites$N_sites[match(unlist(active_genes), gene_sites$gene_id)],
        na.rm = TRUE
      ),
      S_sites = sum(
        gene_sites$S_sites[match(unlist(active_genes), gene_sites$gene_id)],
        na.rm = TRUE
      ),
      dN = (N_count + 0.5) / N_sites,
      dS = (S_count + 0.5) / S_sites,
      dNdS = dN / dS
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, N_count, S_count, N_sites, S_sites, dNdS)
  
  # ---------- 2) bootstrap：对 dNdS vs 1 做双侧检验 ----------
  set.seed(seed)
  
  sample_boot_sub <- purrr::map_dfr(samples_all, function(smp) {
    
    df_one <- mut_sub %>%
      dplyr::filter(sample == smp) %>%
      dplyr::select(gene_id, is_nonsyn, is_syn)
    
    active_genes <- unique(df_one$gene_id[!is.na(df_one$gene_id)])
    
    N_sites_s <- sum(gene_sites$N_sites[match(active_genes, gene_sites$gene_id)], na.rm = TRUE)
    S_sites_s <- sum(gene_sites$S_sites[match(active_genes, gene_sites$gene_id)], na.rm = TRUE)
    
    n_events <- nrow(df_one)
    if (n_events == 0 || !is.finite(N_sites_s) || !is.finite(S_sites_s) ||
        N_sites_s <= 0 || S_sites_s <= 0) {
      return(tibble::tibble(
        sample = smp,
        p_value = NA_real_,
        p_adj = NA_real_,
        neg_log10_p_adj = NA_real_,
        dNdS_boot_mean = NA_real_
      ))
    }
    
    idx_list <- replicate(
      bootstrap_reps,
      sample.int(n_events, size = n_events, replace = TRUE),
      simplify = FALSE
    )
    
    # 关键修复：避免 purrr::vapply（你当前 purrr 版本没有导出 vapply）
    dNdS_boot_vec <- sapply(idx_list, function(idx) {
      boot_dat <- df_one[idx, , drop = FALSE]
      N_count_boot <- sum(boot_dat$is_nonsyn, na.rm = TRUE)
      S_count_boot <- sum(boot_dat$is_syn,    na.rm = TRUE)
      dN_boot <- (N_count_boot + 0.5) / N_sites_s
      dS_boot <- (S_count_boot + 0.5) / S_sites_s
      dN_boot / dS_boot
    })
    
    dNdS_boot_vec <- dNdS_boot_vec[is.finite(dNdS_boot_vec)]
    if (length(dNdS_boot_vec) == 0) {
      return(tibble::tibble(
        sample = smp,
        p_value = NA_real_,
        p_adj = NA_real_,
        neg_log10_p_adj = NA_real_,
        dNdS_boot_mean = NA_real_
      ))
    }
    
    # 双侧 p
    p_upper <- mean(dNdS_boot_vec >= 1, na.rm = TRUE)
    p_lower <- mean(dNdS_boot_vec <= 1, na.rm = TRUE)
    p_value <- 2 * pmin(p_upper, p_lower)
    p_value <- min(p_value, 1)
    
    dNdS_boot_mean <- mean(dNdS_boot_vec, na.rm = TRUE)
    
    tibble::tibble(
      sample = smp,
      p_value = p_value,
      dNdS_boot_mean = dNdS_boot_mean
    )
  })
  
  sample_boot_sub <- sample_boot_sub %>%
    dplyr::mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      p_adj = dplyr::if_else(is.na(p_adj), NA_real_, pmax(p_adj, 1e-300)),
      neg_log10_p_adj = -log10(p_adj)
    )
  
  # ---------- 3) 合并点估计 ----------
  sample_plot_sub <- sample_boot_sub %>%
    dplyr::left_join(sample_point_sub, by = "sample") %>%
    dplyr::mutate(
      signif_class = dplyr::case_when(
        !is.na(p_adj) & p_adj < 0.05 & dNdS > 1 ~ "Significant (dNdS > 1)",
        !is.na(p_adj) & p_adj < 0.05 & dNdS < 1 ~ "Significant (dNdS < 1)",
        TRUE ~ "Not significant"
      )
    ) %>%
    dplyr::select(
      sample,
      N_count, S_count, N_sites, S_sites,
      dNdS,
      dNdS_boot_mean,
      p_value, p_adj, neg_log10_p_adj,
      signif_class
    )
  
  # ---------- 4) 输出 TSV ----------
  label_safe <- dplyr::case_when(
    signif_label == "Significant (dNdS > 1)" ~ "Significant_dNdS_morethan1_",
    signif_label == "Significant (dNdS < 1)" ~ "Significant_dNdS_lessthan1_",
    signif_label == "Not significant"         ~ "Not_significant_",
    TRUE                                       ~ paste0("signif_class_", seq_len(1))
  )
  
  out_tsv <- sprintf(
    "sample_dNdS_vs1_padj_double_sided_by_signif_class_%s.tsv",
    label_safe
  )
  
  readr::write_tsv(sample_plot_sub, out_tsv)
  cat(sprintf("[OK] %s -> %s\n", signif_label, out_tsv))
  
  # ---------- 5) 火山图 ----------
  colors_sample <- c(
    "Significant (dNdS > 1)" = "#E41A1C",
    "Significant (dNdS < 1)" = "#377EB8",
    "Not significant" = "grey60"
  )
  
  p_sub <- ggplot(sample_plot_sub, aes(x = dNdS, y = neg_log10_p_adj, color = signif_class)) +
    geom_point(alpha = 0.85, size = 2.5, na.rm = TRUE) +
    geom_vline(xintercept = 1, linetype = "solid", color = "grey50", linewidth = 0.4) +
    theme_bw(base_size = 12) +
    scale_color_manual(values = colors_sample) +
    labs(
      x = "dNdS (point estimate; cutoff = 1)",
      y = "-log10(p_adj, two-sided; BH corrected)",
      color = "Sample test",
      title = sprintf("Sample-level: dNdS vs 1 (double-sided) | Genes group: %s", signif_label)
    ) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  
  out_png <- sprintf(
    "sample_volcano_dNdS_vs1_padj_double_sided_by_signif_class_%s.png",
    label_safe
  )
  ggsave(out_png, p_sub, width = 10, height = 7, dpi = 300)
  
  cat(sprintf("[OK] %s -> %s\n\n", signif_label, out_png))
  
  sample_plot_sub
}

# ---------- 执行三组 ----------
results_by_group <- list()

for (sg in signif_levels) {
  
  if (!("gene_id" %in% names(gene_dnds_burden))) {
    stop("gene_dnds_burden 中不存在 gene_id 列，请检查列名。")
  }
  if (!("signif_class" %in% names(gene_dnds_burden))) {
    stop("gene_dnds_burden 中不存在 signif_class 列。")
  }
  
  genes_vec <- gene_dnds_burden %>%
    dplyr::filter(signif_class == sg) %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  # 与 gene_sites 取交集（稳健）
  genes_vec <- intersect(genes_vec, gene_sites$gene_id)
  
  if (length(genes_vec) == 0) {
    warning(sprintf("signif_class=%s：genes_vec 与 gene_sites 无交集，跳过。", sg))
    next
  }
  
  results_by_group[[sg]] <- run_sample_dnds_vs1_for_gene_set(genes_vec, sg)
}

cat("\n[OK] 三组 signif_class 的样本级分组分析完成。\n")













# ============================================================
# 20) 基于 PM2.5 及其成分 High/Low 分组的基因 Di 差异分析及组合火山图 (EnhancedVolcano版)
# ============================================================

library(readxl)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(patchwork) # 用于图形拼接
library(EnhancedVolcano) # 用于绘制高级火山图

# 1. 读取包含分组变量的临床/暴露数据
cat("\n[进度] 开始读取暴露分组数据...\n")
exposure_data_path <- "/home/data/t040512/UnfixedMutation/phylogene_comprehensive.xlsx"
exposure_data <- read_excel(exposure_data_path)

# 统一样本列名
sample_col_name <- grep("^Sample_ID$|^Sample$|^sample$", colnames(exposure_data), ignore.case = TRUE)
if (length(sample_col_name) == 0) {
  stop("未在 exposure_data 中找到 Sample_ID 或 Sample 列。")
}
colnames(exposure_data)[sample_col_name[1]] <- "sample"
exposure_data$sample <- as.character(exposure_data$sample)

# 2. 提取所有需要分析的分组变量名 (144个)
group_vars <- grep("_group$", colnames(exposure_data), value = TRUE)
cat(sprintf("[进度] 共检测到 %d 个暴露分组变量，准备进行循环分析。\n", length(group_vars)))

# 3. 将突变数据与暴露分组合并
mut_sample_col <- grep("^Sample_ID$|^Sample$|^sample$", colnames(mut_coding), ignore.case = TRUE)
if (length(mut_sample_col) > 0) {
  colnames(mut_coding)[mut_sample_col[1]] <- "sample"
  mut_coding$sample <- as.character(mut_coding$sample)
}

mut_coding_merged <- mut_coding %>%
  left_join(exposure_data %>% select(sample, all_of(group_vars)), by = "sample") %>%
  filter(!is.na(gene_id))

# 4. 创建输出文件夹
output_dir <- "Di_Diff_Analysis_by_Exposure_Groups"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 5. 定义分析函数
analyze_di_diff <- function(df, group_var, gene_sites_df) {
  df_sub <- df %>% filter(!is.na(.data[[group_var]]))
  if (nrow(df_sub) == 0) return(NULL)
  
  gene_counts <- df_sub %>%
    group_by(gene_id, .data[[group_var]]) %>%
    summarise(Mi = n(), .groups = "drop") %>%
    pivot_wider(names_from = .data[[group_var]], values_from = Mi, values_fill = 0)
  
  if (!"High" %in% colnames(gene_counts)) gene_counts$High <- 0
  if (!"Low" %in% colnames(gene_counts)) gene_counts$Low <- 0
  
  N_high <- sum(gene_counts$High, na.rm = TRUE)
  N_low <- sum(gene_counts$Low, na.rm = TRUE)
  
  if (N_high == 0 || N_low == 0) return(NULL)
  
  res <- gene_sites_df %>%
    left_join(gene_counts, by = "gene_id") %>%
    mutate(
      High = replace_na(High, 0),
      Low = replace_na(Low, 0)
    ) %>%
    filter(!is.na(Li_bp), Li_bp > 0)
  
  res <- res %>%
    mutate(
      # 严格按照参考文献：Di = Mi / (N * Li)
      Di_High = High / (N_high * Li_bp),
      Di_Low = Low / (N_low * Li_bp),
      # 使用 0.5 伪计数计算相对 Fold Change，避免极小值精度问题
      # FC_Di = (Di_High / Di_Low) = (High / N_high) / (Low / N_low) = (High * N_low) / (Low * N_high)
      FC_Di = ((High + 0.5) * N_low) / ((Low + 0.5) * N_high),
      log2FC_Di = log2(FC_Di)
    )
  
  # 处理可能的 Inf 值
  res$log2FC_Di[is.infinite(res$log2FC_Di)] <- NA
  
  res <- res %>%
    mutate(
      p_value = pmap_dbl(list(High, Low), function(h, l) {
        mat <- matrix(c(h, N_high - h, l, N_low - l), nrow = 2, byrow = FALSE)
        tryCatch(fisher.test(mat)$p.value, error = function(e) NA_real_)
      }),
      OR = pmap_dbl(list(High, Low), function(h, l) {
        mat <- matrix(c(h, N_high - h, l, N_low - l), nrow = 2, byrow = FALSE)
        tryCatch(as.numeric(fisher.test(mat)$estimate), error = function(e) NA_real_)
      })
    )
  
  res$padj <- p.adjust(res$p_value, method = "BH")
  res$neg_log10_padj <- -log10(res$padj + 1e-300)
  
  res <- res %>%
    mutate(
      signif = case_when(
        !is.na(padj) & padj < 0.05 & log2FC_Di > 0 ~ "Up in High",
        !is.na(padj) & padj < 0.05 & log2FC_Di < 0 ~ "Down in High",
        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(padj)
  
  return(res)
}

# 6. 循环遍历所有 144 个变量，执行分析、画图并保存
all_results_list <- list()
plot_list <- list() # 用于存储所有单个火山图

# 定义污染物列表，用于最后拼接
pollutants <- c("PM2.5", "SO4", "NO3", "NH4", "OM", "BC")

cat("\n[进度] 开始批量计算与绘图...\n")
for (gvar in group_vars) {
  
  cat(sprintf("  正在处理: %s ...\n", gvar))
  res_df <- analyze_di_diff(mut_coding_merged, gvar, gene_sites)
  if (is.null(res_df)) next
  
  all_results_list[[gvar]] <- res_df
  
  # 提取暴露期长度作为子图标题 (例如 w1_30)
  window_name <- gvar
  for(p in pollutants) {
    window_name <- sub(paste0("^", p, "_"), "", window_name)
  }
  window_name <- gsub("_group$", "", window_name)
  
  # 清理含有 NA 的行，防止 EnhancedVolcano 报错
  res_df_clean <- res_df %>% filter(!is.na(log2FC_Di), !is.na(padj))
  
  # 【关键修复1】：将空的基因标签替换为 "Unannotated"，防止画空白标签
  res_df_clean$gene_id[is.na(res_df_clean$gene_id) | res_df_clean$gene_id == ""] <- "Unannotated"
  
  # 如果清理后没有数据，跳过绘图
  if (nrow(res_df_clean) == 0) {
    plot_list[[gvar]] <- ggplot() + theme_void() + ggtitle(paste0(window_name, "\n(No valid data)"))
    next
  }
  
  # 绘制单个 EnhancedVolcano 火山图 (按照参考代码风格配置)
  p <- EnhancedVolcano(res_df_clean,
                       lab = res_df_clean$gene_id,
                       x = 'log2FC_Di',
                       y = 'padj',
                       xlab = bquote(~Log[2]~ 'fold change (Di)'),
                       ylab = bquote(~-log[10]~ '(FDR)'),
                       pCutoff = 0.05,       # 保留 0.05 适合 FDR
                       FCcutoff = 1.0,       # 【关键修复2】：设为 1.0 (即2倍变化)，让更多基因显示标签
                       cutoffLineType = 'twodash', # 添加虚线样式
                       cutoffLineWidth = 0.8,      # 添加线宽
                       pointSize = 2.0,      # 适当缩小以适应组合图
                       labSize = 3.5,        # 适当缩小以适应组合图
                       colAlpha = 1,         # 按照参考代码设为不透明
                       legendLabels = c('Not sig.','Log (base 2) FC','p-value','p-value & Log (base 2) FC'),
                       legendPosition = 'right', # 按照参考代码显示图例
                       legendLabSize = 10,   # 适当缩小图例字体
                       legendIconSize = 3.0, # 适当缩小图例图标
                       # 添加连接线配置
                       drawConnectors = TRUE,
                       widthConnectors = 0.75, # 按照参考代码设为 0.75
                       max.overlaps = 50,    # 【关键修复3】：增加最大重叠数，防止 ggrepel 隐藏标签
                       title = window_name
  ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.text = element_text(size = 7),
      plot.margin = margin(5, 5, 5, 5) # 适当增加边距给连接线留空间
    )
  
  plot_list[[gvar]] <- p
  
  # 提取显著差异基因单独存为一个 tsv
  sig_genes <- res_df %>% filter(padj < 0.05)
  if (nrow(sig_genes) > 0) {
    sig_tsv <- file.path(output_dir, paste0("Significant_Genes_", gvar, ".tsv"))
    write_tsv(sig_genes, sig_tsv)
  }
}

# 7. 拼接并保存组合火山图 (每种污染物一张图，4列6行，严格按时间排序)
cat("\n[进度] 正在拼接各污染物的组合火山图...\n")
for (p_name in pollutants) {
  pattern <- paste0("^", p_name, "_")
  current_plots <- plot_list[grepl(pattern, names(plot_list))]
  
  if (length(current_plots) > 0) {
    # 【关键修复4】：修正正则表达式，兼容 w1_ 和 window1_，严格提取末尾数字排序
    window_nums <- as.numeric(gsub(".*_(\\d+)_group$", "\\1", names(current_plots)))
    current_plots <- current_plots[order(window_nums)]
    
    combined_plot <- wrap_plots(current_plots, ncol = 4) +
      plot_annotation(
        title = paste0("Gene Di Differential Analysis: ", p_name),
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
    
    # 【关键修复5】：画布尺寸放大至 16x20，以容纳连接线和标签
    combined_filename <- file.path(output_dir, paste0("Combined_Volcano_", p_name, "_Enhanced.tiff"))
    ggsave(combined_filename, combined_plot, width = 16, height = 20, dpi = 300, compression = "lzw")
    cat(sprintf("  -> %s 组合火山图已保存\n", p_name))
  }
}

# 8. 将所有结果汇总写入一个多页 Excel
cat("\n[进度] 正在汇总输出 Excel 结果表...\n")
final_excel_path <- file.path(output_dir, "All_Genes_Di_Diff_Results.xlsx")
wb <- createWorkbook()

for (gvar in names(all_results_list)) {
  sheet_name <- gsub("_window1_", "_w1_", gvar)
  sheet_name <- substr(sheet_name, 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, all_results_list[[gvar]])
}
saveWorkbook(wb, final_excel_path, overwrite = TRUE)
cat(sprintf("汇总 Excel 表已保存至: %s\n", final_excel_path))

cat("\n全部脚本运行结束。\n")











# ============================================================
# 21) 显著差异基因的功能注释 (直接从 H37Rv.gb 提取注释)
# ============================================================

library(readxl)
library(openxlsx)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(readr)

# 1. 设置文件路径
setwd("/home/data/t040512/UnfixedMutation")
gb_file <- "H37Rv.gb"
input_excel_path <- "Di_Diff_Analysis_by_Exposure_Groups/All_Genes_Di_Diff_Results.xlsx"
output_excel_path <- "Di_Diff_Analysis_by_Exposure_Groups/Significant_Genes_Annotated.xlsx"

# 2. 解析 H37Rv.gb 文件，提取基因注释信息
cat("[进度] 正在解析 H37Rv.gb 提取基因注释...\n")
gb_lines <- read_lines(gb_file)
feat_start <- which(str_detect(gb_lines, "^FEATURES"))
origin_start <- which(str_detect(gb_lines, "^ORIGIN"))
feat_block <- gb_lines[(feat_start + 1):(origin_start - 1)]
is_feat_head <- str_detect(feat_block, "^ {5}\\S")
head_idx <- which(is_feat_head)

parse_annotation <- function(start_i, end_i) {
  block <- feat_block[start_i:end_i]
  ftype <- str_trim(str_sub(block[1], 6, 21))
  if (ftype != "CDS") return(NULL)
  joined <- paste(block, collapse = " ")
  locus_tag <- str_match(joined, '/locus_tag="([^"]+)"')[, 2]
  gene_name <- str_match(joined, '/gene="([^"]+)"')[, 2]
  product <- str_match(joined, '/product="([^"]+)"')[, 2]
  if (is.na(locus_tag)) return(NULL)
  tibble(gene_id = locus_tag, Gene_Name = ifelse(is.na(gene_name), "", gene_name), Product = ifelse(is.na(product), "", product))
}

gene_annotation <- map2_dfr(head_idx, c(head_idx[-1] - 1, length(feat_block)), parse_annotation) %>% distinct(gene_id, .keep_all = TRUE)
cat(sprintf("[进度] 成功提取 %d 个基因的注释信息。\n", nrow(gene_annotation)))

# 3. 提取显著基因并合并注释
cat("[进度] 正在读取差异分析结果并进行注释...\n")
sheet_names <- getSheetNames(input_excel_path)
annotated_sig_list <- list()

for (sheet in sheet_names) {
  res_df <- readWorkbook(input_excel_path, sheet = sheet)
  sig_df <- res_df %>% filter(padj < 0.05)
  if (nrow(sig_df) > 0) {
    sig_annotated <- sig_df %>%
      left_join(gene_annotation, by = "gene_id") %>%
      mutate(Display_Name = ifelse(Gene_Name == "" | is.na(Gene_Name), gene_id, Gene_Name)) %>%
      select(gene_id, Gene_Name, Display_Name, Product, signif, log2FC_Di, padj, everything())
    annotated_sig_list[[sheet]] <- sig_annotated
  }
}

# 4. 输出带有注释的 Excel 文件
if (length(annotated_sig_list) > 0) {
  wb_sig <- createWorkbook()
  header_style <- createStyle(fontSize = 11, fontColour = "#FFFFFF", fgFill = "#4F81BD", halign = "center", textDecoration = "bold")
  
  for (sheet in names(annotated_sig_list)) {
    sheet_name <- substr(sheet, 1, 31)
    addWorksheet(wb_sig, sheet_name)
    writeData(wb_sig, sheet_name, annotated_sig_list[[sheet]])
    addStyle(wb_sig, sheet_name, header_style, rows = 1, cols = 1:ncol(annotated_sig_list[[sheet]]), gridExpand = TRUE)
    freezePane(wb_sig, sheet_name, firstRow = TRUE)
    setColWidths(wb_sig, sheet_name, cols = 1:ncol(annotated_sig_list[[sheet]]), widths = "auto")
  }
  
  # 全局统计
  all_sig_combined <- bind_rows(annotated_sig_list, .id = "Exposure_Group")
  summary_stats <- all_sig_combined %>% group_by(Exposure_Group, signif) %>% summarise(Count = n(), .groups = "drop") %>% arrange(Exposure_Group, signif)
  
  addWorksheet(wb_sig, "Summary_Stats", gridLines = FALSE)
  writeData(wb_sig, "Summary_Stats", summary_stats)
  addStyle(wb_sig, "Summary_Stats", header_style, rows = 1, cols = 1:ncol(summary_stats), gridExpand = TRUE)
  freezePane(wb_sig, "Summary_Stats", firstRow = TRUE)
  worksheetOrder(wb_sig) <- c(which(names(wb_sig) == "Summary_Stats"), setdiff(seq_along(names(wb_sig)), which(names(wb_sig) == "Summary_Stats")))
  
  saveWorkbook(wb_sig, output_excel_path, overwrite = TRUE)
  cat(sprintf("\n[完成] 显著差异基因的功能注释已保存至: %s\n", output_excel_path))
} else {
  cat("\n[提示] 在所有暴露窗口中均未发现 padj < 0.05 的显著差异基因。\n")
}












# ============================================================
# 22) 绘制各暴露窗口显著差异基因的柱形图及 4x6 组合大图 (倾斜45度)
# ============================================================

library(ggplot2)
library(ggsci)    # 提供 jco 期刊配色
library(openxlsx)
library(dplyr)
library(stringr)
library(patchwork)

# 1. 设置路径与读取数据
setwd("/home/data/t040512/UnfixedMutation")
input_excel_path <- "Di_Diff_Analysis_by_Exposure_Groups/Significant_Genes_Annotated.xlsx"
output_dir <- "Di_Diff_Analysis_by_Exposure_Groups/Barplots_Sig_Genes"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. 读取包含显著基因注释的 Excel
cat("[进度] 正在读取显著差异基因数据...\n")
sheet_names <- getSheetNames(input_excel_path)
sheet_names <- sheet_names[sheet_names != "Summary_Stats"]

plot_list <- list()

# 3. 循环遍历每个暴露组，绘制单个柱形图
for (sheet in sheet_names) {
  
  df <- readWorkbook(input_excel_path, sheet = sheet)
  
  if (nrow(df) == 0) next
  
  df <- df %>%
    mutate(
      Display_Name = ifelse(is.na(Display_Name) | Display_Name == "", as.character(gene_id), as.character(Display_Name)),
      Direction = ifelse(log2FC_Di > 0, "Up in High", "Down in High"),
      # 按照 log2FC_Di 升序排列基因
      Display_Name = factor(Display_Name, levels = unique(Display_Name[order(df$log2FC_Di)]))
    )
  
  # 使用原生 ggplot2 绘制柱形图
  p <- ggplot(df, aes(x = Display_Name, y = log2FC_Di, fill = Direction)) +
    geom_col(color = "white", linewidth = 0.2) +  
    scale_fill_jco() +                             
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    labs(
      title = sheet,
      y = "log2 Fold Change (Di)",
      x = NULL,
      fill = "Direction"
    ) +
    scale_y_continuous(expand = c(0, 0)) +         
    theme_classic() +                              
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      # 倾斜45度并右对齐
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7), 
      axis.text.y = element_text(size = 8),
      legend.position = "top"                      
    )
  
  plot_list[[sheet]] <- p
  
  # 保存单个柱形图
  plot_width <- max(5, nrow(df) * 0.25 + 2)
  ggsave(filename = file.path(output_dir, paste0("Barplot_", sheet, ".tiff")), 
         plot = p, width = plot_width, height = 4, dpi = 300, compression = "lzw")
}

# 4. 按污染物拼接组合柱形图 (严格 4列6行 画布布局)
cat("\n[进度] 正在拼接各污染物的 4x6 组合柱形图...\n")
pollutants <- c("PM2.5", "SO4", "NO3", "NH4", "OM", "BC")

# 动态获取所有窗口名并去重
all_exposure_windows <- unique(sub("^(PM2.5|SO4|NO3|NH4|OM|BC)_", "", names(plot_list)))

for (p_name in pollutants) {
  current_plots <- list()
  
  # 遍历所有暴露窗口，确保按顺序填充 4x6=24 个位置
  for (window in all_exposure_windows) {
    sheet_name <- paste0(p_name, "_", window)
    
    if (sheet_name %in% names(plot_list)) {
      current_plots[[window]] <- plot_list[[sheet_name]]
    } else {
      # 空白占位图
      current_plots[[window]] <- ggplot() + theme_void() + 
        labs(title = sheet_name) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10, color = "grey50"))
    }
  }
  
  if (length(current_plots) > 0) {
    # 修改此处：强制 ncol=4, nrow=6
    combined_plot <- wrap_plots(current_plots, ncol = 4, nrow = 6) +
      plot_annotation(
        title = paste0("Significant Gene Di Changes: ", p_name),
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
    
    combined_filename <- file.path(output_dir, paste0("Combined_Barplot_", p_name, "_4x6.tiff"))
    # 修改此处：调整画布大小，更宽更扁 (宽16英寸，高18英寸)
    ggsave(combined_filename, combined_plot, width = 16, height = 18, dpi = 300, compression = "lzw")
    cat(sprintf("  -> %s 4x6 组合柱形图已保存\n", p_name))
  }
}

cat("\n[完成] 所有柱形图已生成并保存至文件夹: ", output_dir, "\n")









# ============================================================
# 23) 基于 MycoBrowser 功能分类的富集分析及可视化 (Release 5 完整版)
# ============================================================

library(openxlsx)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)

# 1. 设置路径与读取数据
setwd("/home/data/t040512/UnfixedMutation")
input_excel_path <- "Di_Diff_Analysis_by_Exposure_Groups/All_Genes_Di_Diff_Results.xlsx"
output_dir <- "Di_Diff_Analysis_by_Exposure_Groups/MycoBrowser_Enrichment"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. 读取并处理 MycoBrowser 注释文件 (针对 Release 5 真实格式)
myco_file <- "Mycobacterium_tuberculosis_H37Rv_txt_v5.txt"
if (!file.exists(myco_file)) {
  stop("未找到 MycoBrowser 注释文件！请检查文件是否在工作目录下。")
}

cat("[进度] 正在读取 MycoBrowser Release 5 注释文件并构建二级分类...\n")

# 读取 TSV 文件，并清洗列名
myco_anno <- read_tsv(myco_file, show_col_types = FALSE) %>%
  rename_with(~ tolower(gsub("[^a-zA-Z0-9_]", "_", .x))) %>%
  select(gene_id = locus, gene_name = name, function_desc = `function`, product, functional_category) %>%
  filter(!is.na(functional_category), functional_category != "")

# 清洗空格并统一转为大写，确保与差异分析结果的基因ID格式完全一致
myco_anno$gene_id <- str_trim(myco_anno$gene_id) %>% toupper()
myco_anno$functional_category <- str_trim(myco_anno$functional_category) 
myco_anno$gene_name <- ifelse(myco_anno$gene_name == "" | is.na(myco_anno$gene_name), myco_anno$gene_id, myco_anno$gene_name)

# --- 核心新增：基于 Function 和 Product 描述匹配二级精细分类 ---
myco_anno <- myco_anno %>%
  mutate(
    func_lower = tolower(paste(function_desc, product, sep = " ")),
    detailed_category = case_when(
      str_detect(func_lower, "dna|replication|repair|recombination|rna polymerase|helicase|ligase|gyrase|nucleotide excision|mismatch|transcription") ~ "DNA replication, repair, recombination and transcription",
      str_detect(func_lower, "regulator|repressor|activator|two-component|sigma factor|transcriptional reg") ~ "Regulatory proteins",
      str_detect(func_lower, "cell wall|membrane|transport|porin|lps|murein|peptidoglycan|divis|secretion") ~ "Cell wall and cell processes",
      str_detect(func_lower, "metabolism|dehydrogenase|oxidoreductase|tca|respiration|atp synthase|kinase|isomerase|hydrolase|transferase") ~ "Intermediary metabolism and respiration",
      str_detect(func_lower, "virulence|toxin|antitoxin|stress|heat shock|chaperone|detoxif") ~ "Virulence, detoxification, adaptation",
      TRUE ~ functional_category
    )
  )

cat(sprintf("[进度] 成功加载 %d 个基因的 MycoBrowser 功能注释（含二级分类）。\n", nrow(myco_anno)))

# 3. 定义 Fisher 富集分析函数
perform_myco_enrichment <- function(sig_genes, bg_genes, anno_df, category_col) {
  sig_anno <- anno_df %>% filter(gene_id %in% sig_genes)
  bg_anno <- anno_df %>% filter(gene_id %in% bg_genes)
  
  if (nrow(sig_anno) == 0) return(NULL)
  
  categories <- unique(anno_df[[category_col]])
  
  res_list <- lapply(categories, function(cat) {
    bg_in_cat <- sum(bg_anno[[category_col]] == cat)
    bg_total <- nrow(bg_anno)
    sig_in_cat <- sum(sig_anno[[category_col]] == cat)
    sig_total <- nrow(sig_anno)
    
    if (bg_in_cat == 0) return(NULL)
    
    mat <- matrix(c(sig_in_cat, sig_total - sig_in_cat,
                    bg_in_cat, bg_total - bg_in_cat),
                  nrow = 2, byrow = TRUE)
    
    ft <- tryCatch(fisher.test(mat, alternative = "greater"), error = function(e) NULL)
    if (is.null(ft)) return(NULL)
    
    data.frame(
      Category_Level = ifelse(category_col == "functional_category", "Level1_Broad", "Level2_Detailed"),
      Category = cat,
      Sig_in_cat = sig_in_cat,
      Sig_total = sig_total,
      Bg_in_cat = bg_in_cat,
      Bg_total = bg_total,
      FoldEnrichment = (sig_in_cat / sig_total) / (bg_in_cat / bg_total),
      pvalue = ft$p.value
    )
  })
  
  res_df <- bind_rows(res_list) %>% filter(Sig_in_cat > 0)
  return(res_df)
}

# 4. 遍历所有暴露组进行富集分析
cat("[进度] 开始基于 MycoBrowser 的功能富集分析...\n")
sheet_names <- getSheetNames(input_excel_path)
all_myco_results <- list()

# 获取背景基因并统一格式
bg_genes <- unique(unlist(lapply(sheet_names, function(sn) readWorkbook(input_excel_path, sheet = sn)$gene_id)))
bg_genes <- str_trim(bg_genes) %>% toupper()

for (sheet in sheet_names) {
  res_df <- readWorkbook(input_excel_path, sheet = sheet)
  # 统一差异基因格式：去空格、转大写
  sig_genes <- res_df %>% filter(padj < 0.05) %>% pull(gene_id) %>% str_trim() %>% toupper() %>% unique()
  
  cat(sprintf("  -> 处理组: %s | 显著基因数: %d\n", sheet, length(sig_genes)))
  
  # 修改点1：放宽限制，只要有1个显著基因就尝试分析，避免低突变组被遗漏
  if (length(sig_genes) < 1) next
  
  enrich_res_l1 <- perform_myco_enrichment(sig_genes, bg_genes, myco_anno, "functional_category")
  enrich_res_l2 <- perform_myco_enrichment(sig_genes, bg_genes, myco_anno, "detailed_category")
  
  enrich_res <- bind_rows(enrich_res_l1, enrich_res_l2)
  
  if (!is.null(enrich_res) && nrow(enrich_res) > 0) {
    enrich_res$Exposure_Group <- sheet
    all_myco_results[[sheet]] <- enrich_res
    cat(sprintf("     成功富集到 %d 条结果 (L1: %d, L2: %d)\n", nrow(enrich_res), nrow(enrich_res_l1), nrow(enrich_res_l2)))
  } else {
    cat(sprintf("     未富集到任何结果。\n"))
  }
}

# 5. 汇总结果、输出 Excel 与绘图
if (length(all_myco_results) > 0) {
  final_myco_df <- bind_rows(all_myco_results)
  
  final_myco_df <- final_myco_df %>%
    mutate(
      Pollutant = str_extract(Exposure_Group, "^(PM2\\.5|SO4|NO3|NH4|OM|BC)"),
      neg_log10_pvalue = -log10(pvalue + 1e-300)
    )
  
  # 保存 Excel
  myco_excel_path <- file.path(output_dir, "MycoBrowser_Enrichment_Results.xlsx")
  wb_myco <- createWorkbook()
  addWorksheet(wb_myco, "MycoBrowser_Enrichment")
  writeData(wb_myco, "MycoBrowser_Enrichment", final_myco_df)
  
  header_style <- createStyle(fontSize = 11, fontColour = "#FFFFFF", fgFill = "#4F81BD", halign = "center", textDecoration = "bold")
  addStyle(wb_myco, "MycoBrowser_Enrichment", header_style, rows = 1, cols = 1:ncol(final_myco_df), gridExpand = TRUE)
  freezePane(wb_myco, "MycoBrowser_Enrichment", firstRow = TRUE)
  saveWorkbook(wb_myco, myco_excel_path, overwrite = TRUE)
  cat(sprintf("[完成] MycoBrowser 富集分析结果已保存至: %s\n", myco_excel_path))
  
  # 绘制柱形图函数
  draw_myco_bar <- function(df, group_name) {
    # 修改点2：优先画二级分类，如果没有则回退画一级分类，防止图空白
    df_plot <- df %>% filter(Category_Level == "Level2_Detailed")
    
    if (nrow(df_plot) == 0) {
      df_plot <- df %>% filter(Category_Level == "Level1_Broad")
    }
    
    if (nrow(df_plot) == 0) return(NULL)
    
    df_plot <- df_plot %>% 
      arrange(FoldEnrichment) %>% 
      mutate(Category = factor(Category, levels = Category))
    
    text_size <- ifelse(nrow(df_plot) > 10, 6, 8)
    
    p <- ggplot(df_plot, aes(x = FoldEnrichment, y = Category, fill = neg_log10_pvalue)) +
      geom_col(width = 0.7, color = "black", linewidth = 0.1) +
      geom_text(aes(label = Sig_in_cat), hjust = -0.1, size = text_size * 0.35, color = "grey30") +
      scale_fill_gradient(low = "#B3CDE3", high = "#D95F0E") +
      theme_bw(base_size = text_size) +
      labs(x = "Fold Enrichment", y = NULL, title = group_name, fill = "-log10(P)") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = text_size + 1),
        axis.text.y = element_text(color = "black", face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.5),
        legend.position = "right", plot.margin = margin(5, 25, 5, 5)
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.3)))
    return(p)
  }
  
  plot_list <- list()
  groups_with_results <- unique(final_myco_df$Exposure_Group)
  for (grp in groups_with_results) {
    df_sub <- final_myco_df %>% filter(Exposure_Group == grp)
    p <- draw_myco_bar(df_sub, grp)
    if (!is.null(p)) plot_list[[grp]] <- p
  }
  
  time_windows <- c(30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 
                    390, 420, 450, 480, 510, 540, 570, 600, 630, 660, 690, 720)
  
  pollutants <- c("PM2.5", "SO4", "NO3", "NH4", "OM", "BC")
  cat("[进度] 正在拼接各污染物的 MycoBrowser 组合柱形图 (3x8严格排序布局)...\n")
  
  for (p_name in pollutants) {
    current_plots <- list()
    
    for (tw in time_windows) {
      grp_name <- paste0(p_name, "_w1_", tw, "_group")
      
      if (grp_name %in% names(plot_list)) {
        current_plots[[grp_name]] <- plot_list[[grp_name]]
      } else {
        # 修改点3：如果该组确实无富集结果，生成带提示的空白占位图
        current_plots[[grp_name]] <- ggplot() + theme_void() + 
          labs(title = grp_name, subtitle = "No significant enrichment") + 
          theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", color = "grey50"),
                plot.subtitle = element_text(hjust = 0.5, size = 6, color = "red"))
      }
    }
    
    combined_plot <- wrap_plots(current_plots, ncol = 3, nrow = 8) +
      plot_annotation(title = paste0("MycoBrowser Detailed Functional Enrichment: ", p_name), 
                      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))
    
    combined_filename <- file.path(output_dir, paste0("Combined_MycoBrowser_Enrichment_", p_name, "_3x8.tiff"))
    ggsave(combined_filename, combined_plot, width = 16, height = 20, dpi = 300, compression = "lzw")
    cat(sprintf("  -> %s 3x8 组合柱形图已保存\n", p_name))
  }
  cat("[完成] 所有 MycoBrowser 富集分析柱形图已生成并保存。\n")
} else {
  cat("\n[提示] 所有暴露组均未在 MycoBrowser 分类中富集到显著结果。\n")
}






















# ============================================================
# 24) 绘制 MycoBrowser 富集分析高级组合图 (柳叶刀风格)
# ============================================================

library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(patchwork)
library(readxl)

# 1. 设置工作目录并读取数据 (修正了文件路径)
setwd("/home/data/t040512/UnfixedMutation")
# 明确指定文件在子目录下
file_path <- "Di_Diff_Analysis_by_Exposure_Groups/MycoBrowser_Enrichment/MycoBrowser_Enrichment_Results.xlsx"

if (!file.exists(file_path)) {
  stop("未找到富集分析结果文件，请检查路径！")
}

df_enrich <- read_excel(file_path)

# 2. 筛选显著结果用于绘图
df_sig <- df_enrich %>% filter(pvalue < 0.05)

if (nrow(df_sig) > 0) {
  # 3. 定义柳叶刀风格主题
  lancet_theme <- theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.y = element_text(color = "black", face = "bold"),
      axis.text.x = element_text(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.5),
      legend.position = "right",
      plot.margin = margin(5, 15, 5, 5)
    )
  
  # 4. 绘制分面柱状图函数
  draw_lancet_bar <- function(df, group_name) {
    df_plot <- df %>% 
      arrange(FoldEnrichment) %>% 
      mutate(Category = factor(Category, levels = Category))
    
    text_size <- ifelse(nrow(df_plot) > 10, 6, 8)
    
    p <- ggplot(df_plot, aes(x = FoldEnrichment, y = Category, fill = neg_log10_pvalue)) + # 修改为 neg_log10_pvalue
      geom_col(width = 0.7, color = "black", linewidth = 0.2) +
      # 去掉 "Count: " 前缀，直接显示数字
      geom_text(aes(label = Sig_in_cat), hjust = -0.1, size = text_size * 0.35, color = "grey30") +
      scale_fill_gradient(low = "#B3CDE3", high = "#D95F0E") +
      labs(x = "Fold Enrichment", y = NULL, title = group_name, fill = "-log10(P)") + # 修改为 -log10(P)
      lancet_theme +
      scale_x_continuous(expand = expansion(mult = c(0, 0.3)))
    return(p)
  }
  
  # 5. 循环绘图并拼接
  plot_list <- list()
  groups_with_results <- unique(df_sig$Exposure_Group)
  
  for (grp in groups_with_results) {
    df_sub <- df_sig %>% filter(Exposure_Group == grp)
    if (nrow(df_sub) > 0) { 
      plot_list[[grp]] <- draw_lancet_bar(df_sub, grp) 
    }
  }
  
  pollutants <- c("PM2.5", "SO4", "NO3", "NH4", "OM", "BC")
  cat("\n[进度] 正在拼接各污染物的柳叶刀风格组合图 (3x2布局)...\n")
  
  for (p_name in pollutants) {
    pattern <- paste0("^", p_name, "_")
    current_plots <- plot_list[grepl(pattern, names(plot_list))]
    
    if (length(current_plots) > 0) {
      current_plots <- current_plots[order(names(current_plots))]
      combined_plot <- wrap_plots(current_plots, ncol = 2) +
        plot_annotation(
          title = paste0("MycoBrowser Functional Enrichment: ", p_name), 
          theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
        )
      
      combined_filename <- file.path("Di_Diff_Analysis_by_Exposure_Groups/MycoBrowser_Enrichment", paste0("Fig_MycoBrowser_Enrichment_", p_name, "_Lancet.tiff"))
      ggsave(combined_filename, combined_plot, width = 15, height = 20, dpi = 300, compression = "lzw")
      cat(sprintf("  -> %s 柳叶刀风格组合图已保存\n", p_name))
    }
  }
  cat("\n[完成] 所有柳叶刀风格富集分析组合图已生成。\n")
} else {
  cat("\n[提示] 未找到 pvalue < 0.05 的结果用于绘图。\n")
}







# ============================================================
# 25) 污染物暴露水平与氧化应激碱基突变类型 (C>T, G>A, G>T, C>A) 关系分析
# ============================================================

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(openxlsx)
library(ggplot2)
library(patchwork)
library(janitor)

# 1. 设置输出目录
setwd("/home/data/t040512/UnfixedMutation")
output_dir_ox <- "Oxidative_Mutation_by_Exposure_Groups"
if (!dir.exists(output_dir_ox)) dir.create(output_dir_ox, recursive = TRUE)

# 2. 重新读取原始突变数据，提取 REF 和 ALT 碱基
cat("\n[进度 25.1] 读取原始突变数据并提取碱基改变类型...\n")
mut_raw_ox <- read_excel(mut_xlsx) |> janitor::clean_names()

cat("  [调试] 数据列名: ", paste(names(mut_raw_ox), collapse = ", "), "\n")

# 自动搜索包含碱基改变模式的列 (如 "C>T")
mut_col <- NULL
for (col in names(mut_raw_ox)) {
  if (any(str_detect(as.character(mut_raw_ox[[col]]), "[ACGT]>[ACGT]"))) {
    mut_col <- col
    break
  }
}

if (!is.null(mut_col)) {
  cat(sprintf("  [信息] 在列 '%s' 中检测到碱基改变模式。\n", mut_col))
  mut_raw_ox <- mut_raw_ox %>%
    mutate(mut_pattern = str_extract(str_to_upper(as.character(.data[[mut_col]])), "[ACGT]>[ACGT]"))
} else {
  # 如果没有直接包含 C>T 的列，尝试拼接 ref 和 alt
  ref_col <- grep("^ref$|^ref_base$|^refbase$|^ref_allele$|^reference$", names(mut_raw_ox), ignore.case = TRUE, value = TRUE)[1]
  alt_col <- grep("^alt$|^alt_base$|^altbase$|^alt_allele$|^alternate$|^var$|^variant$", names(mut_raw_ox), ignore.case = TRUE, value = TRUE)[1]
  
  if (!is.na(ref_col) && !is.na(alt_col)) {
    cat(sprintf("  [信息] 通过拼接 '%s' 和 '%s' 生成碱基改变。\n", ref_col, alt_col))
    mut_raw_ox <- mut_raw_ox %>%
      mutate(mut_pattern = paste0(str_to_upper(str_trim(as.character(.data[[ref_col]]))), ">", str_to_upper(str_trim(as.character(.data[[alt_col]])))))
  } else {
    stop("未能找到包含碱基改变的列，请检查数据表列名。")
  }
}

# 过滤掉无法识别碱基改变的突变
mut_ox <- mut_raw_ox %>%
  filter(!is.na(mut_pattern)) %>%
  mutate(
    sample = as.character(sample),
    # 定义目标氧化应激突变类型
    # 胞嘧啶脱氨: C>T, G>A
    # 8-氧鸟嘌呤: G>T, C>A
    is_oxidative = mut_pattern %in% c("C>T", "G>A", "G>T", "C>A")
  )

cat(sprintf("[进度 25.2] 成功识别碱基改变 %d 条，其中氧化应激相关突变 %d 条。\n", 
            nrow(mut_ox), sum(mut_ox$is_oxidative)))

# 3. 合并暴露分组数据 (复用第20部分的 exposure_data 和 group_vars)
# 确保样本列名统一
colnames(mut_ox)[grep("^sample$|^sample_id$", colnames(mut_ox), ignore.case = TRUE)[1]] <- "sample"

mut_ox_merged <- mut_ox %>%
  select(sample, mut_pattern, is_oxidative) %>%
  inner_join(exposure_data %>% select(sample, all_of(group_vars)), by = "sample")

# 4. 定义分析函数：计算 High vs Low 组的氧化突变占比、RR 及 Percentage Change
analyze_ox_ratio <- function(df, group_var) {
  df_sub <- df %>% filter(!is.na(.data[[group_var]]))
  if (nrow(df_sub) == 0) return(NULL)
  
  # 按高/低暴露组汇总总突变和氧化突变数
  summary_df <- df_sub %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      Total_Mut = n(),
      Ox_Mut = sum(is_oxidative, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 确保有 High 和 Low 两组
  if (!"High" %in% summary_df[[group_var]] || !"Low" %in% summary_df[[group_var]]) return(NULL)
  
  high_row <- summary_df %>% filter(.data[[group_var]] == "High")
  low_row <- summary_df %>% filter(.data[[group_var]] == "Low")
  
  # 计算非氧化突变数，用于列联表
  non_ox_high <- high_row$Total_Mut - high_row$Ox_Mut
  non_ox_low <- low_row$Total_Mut - low_row$Ox_Mut
  
  # 计算比例 r (防止总突变为0)
  r_high <- high_row$Ox_Mut / max(high_row$Total_Mut, 1)
  r_low <- low_row$Ox_Mut / max(low_row$Total_Mut, 1)
  
  # 计算 RR = r高 / r低
  RR <- r_high / max(r_low, 1e-6)
  
  # Fisher 精确检验
  mat <- matrix(c(high_row$Ox_Mut, non_ox_high, 
                  low_row$Ox_Mut, non_ox_low), nrow = 2, byrow = TRUE)
  ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
  
  if (is.null(ft)) return(NULL)
  
  # 计算 95% 置信区间
  SE_lnRR <- sqrt((1/max(high_row$Ox_Mut, 0.5)) - (1/max(high_row$Total_Mut, 1)) + 
                    (1/max(low_row$Ox_Mut, 0.5)) - (1/max(low_row$Total_Mut, 1)))
  CI_Lower <- exp(log(RR) - 1.96 * SE_lnRR)
  CI_Upper <- exp(log(RR) + 1.96 * SE_lnRR)
  
  # 计算 Percentage Change 及其 CI
  Percentage_change <- (RR - 1) * 100
  Percentage_change_CI_Lower <- (CI_Lower - 1) * 100
  Percentage_change_CI_Upper <- (CI_Upper - 1) * 100
  
  tibble(
    Exposure_Group = group_var,
    High_Total = high_row$Total_Mut,
    High_Ox = high_row$Ox_Mut,
    High_Ratio = r_high,
    Low_Total = low_row$Total_Mut,
    Low_Ox = low_row$Ox_Mut,
    Low_Ratio = r_low,
    RR = RR,
    RR_CI_Lower = CI_Lower,
    RR_CI_Upper = CI_Upper,
    Percentage_change = Percentage_change,
    Percentage_change_CI_Lower = Percentage_change_CI_Lower,
    Percentage_change_CI_Upper = Percentage_change_CI_Upper,
    P_Value = ft$p.value,
    OR = as.numeric(ft$estimate)
  )
}

# 5. 循环遍历所有 144 个暴露窗口进行计算
cat("[进度 25.3] 正在计算各暴露窗口的氧化突变 RR 及 Percentage Change...\n")
ox_results_list <- lapply(group_vars, function(gvar) {
  res <- analyze_ox_ratio(mut_ox_merged, gvar)
  return(res)
})

ox_results_df <- bind_rows(ox_results_list)

# 增加容错：如果没有任何结果，停止后续绘图
if (nrow(ox_results_df) == 0) {
  cat("[警告] 未计算出任何有效的结果。请检查：\n")
  cat("  1. mut_ox_merged 是否为空（样本名是否匹配）？\n")
  cat("  2. 是否所有组的总突变数都为0？\n")
} else {
  ox_results_df <- ox_results_df %>%
    mutate(
      P_adj = p.adjust(P_Value, method = "BH"),
      neg_log10_padj = -log10(P_adj + 1e-300),
      signif = case_when(
        P_adj < 0.05 & Percentage_change > 0 ~ "Up in High (PC > 0)",
        P_adj < 0.05 & Percentage_change < 0 ~ "Down in High (PC < 0)",
        TRUE ~ "Not significant"
      ),
      # 提取污染物名称和时间窗口用于绘图
      Pollutant = str_extract(Exposure_Group, "^(PM2\\.5|SO4|NO3|NH4|OM|BC)"),
      Time_Window = as.numeric(gsub(".*_(\\d+)_group$", "\\1", Exposure_Group))
    ) %>%
    filter(!is.na(Pollutant), !is.na(Time_Window)) %>%
    arrange(Pollutant, Time_Window)
}

# 6. 输出 Excel 结果
ox_excel_path <- file.path(output_dir_ox, "Oxidative_Mutation_Percentage_Change_Results.xlsx")
wb_ox <- createWorkbook()
addWorksheet(wb_ox, "Ox_PC_Results")
writeData(wb_ox, "Ox_PC_Results", ox_results_df)
header_style <- createStyle(fontSize = 11, fontColour = "#FFFFFF", fgFill = "#4F81BD", halign = "center", textDecoration = "bold")
addStyle(wb_ox, "Ox_PC_Results", header_style, rows = 1, cols = 1:ncol(ox_results_df), gridExpand = TRUE)
freezePane(wb_ox, "Ox_PC_Results", firstRow = TRUE)
saveWorkbook(wb_ox, ox_excel_path, overwrite = TRUE)
cat(sprintf("[完成] Percentage Change 计算结果已保存至: %s\n", ox_excel_path))

# 7. 绘制各污染物的 Percentage Change 随时间窗口变化折线图
cat("[进度 25.4] 正在绘制各污染物的 Percentage Change 时间轨迹图 (PDF)...\n")
pollutants <- c("PM2.5", "SO4", "NO3", "NH4", "OM", "BC")
plot_list_ox <- list()

for (p_name in pollutants) {
  df_plot <- ox_results_df %>% filter(Pollutant == p_name)
  if (nrow(df_plot) == 0) next
  
  # 使用 expression 构建带有正确上下标的标题
  p_title_expr <- switch(p_name,
                         "PM2.5" = expression(paste("Oxidative Mutation % Change: ", PM[2.5])),
                         "SO4" = expression(paste("Oxidative Mutation % Change: ", SO[4]^"2−")),
                         "NO3" = expression(paste("Oxidative Mutation % Change: ", NO[3]^"−")),
                         "NH4" = expression(paste("Oxidative Mutation % Change: ", NH[4]^"+")),
                         "OM" = expression(paste("Oxidative Mutation % Change: OM")),
                         "BC" = expression(paste("Oxidative Mutation % Change: BC")),
                         p_name
  )
  
  p <- ggplot(df_plot, aes(x = Time_Window, y = Percentage_change)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_ribbon(aes(ymin = Percentage_change_CI_Lower, ymax = Percentage_change_CI_Upper), fill = "#3C5488FF", alpha = 0.2) +
    geom_line(color = "#3C5488FF", linewidth = 1) +
    geom_point(aes(color = signif, size = signif != "Not significant")) +
    scale_color_manual(values = c("Up in High (PC > 0)" = "#E41A1C", 
                                  "Down in High (PC < 0)" = "#377EB8", 
                                  "Not significant" = "grey60")) +
    scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1.5)) +
    labs(
      title = p_title_expr,
      subtitle = "Mutations due to oxidative stress (cytosine deamination: C>T, G>A; 8-oxoguanine: G>T, C>A)",
      x = "Exposure Time Window (Days)",
      y = "Percentage Change (%) (High vs Low)",
      color = "Significance"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    guides(size = "none")
  
  plot_list_ox[[p_name]] <- p
  
  # 保存单个图 (PDF 格式)
  ggsave(filename = file.path(output_dir_ox, paste0("Percentage_Change_Trajectory_", p_name, ".pdf")), 
         plot = p, width = 8, height = 5, device = cairo_pdf)
}

# 8. 拼接组合大图 (2x3 布局)
if (length(plot_list_ox) > 0) {
  combined_plot_ox <- wrap_plots(plot_list_ox, ncol = 3) +
    plot_annotation(
      title = "Oxidative Stress Mutation Percentage Change by PM2.5 Components",
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
    )
  
  combined_filename_ox <- file.path(output_dir_ox, "Combined_Oxidative_Percentage_Change_Trajectories.pdf")
  ggsave(combined_filename_ox, combined_plot_ox, width = 18, height = 10, device = cairo_pdf)
  cat(sprintf("[完成] 组合轨迹图已保存至: %s\n", combined_filename_ox))
}

cat("\n[全部完成] 第 25 节：氧化应激碱基突变关系分析执行完毕。\n")