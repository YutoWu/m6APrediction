#' DNA序列编码函数
#'
#' 将DNA字符串转换为编码后的数据框
#'
#' @param dna_strings 字符向量，包含DNA序列
#' @return 编码后的数据框
#' @examples
#' dna_encoding("ATCG")
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' 多样本m6A位点预测
#'
#' 该函数使用训练好的机器学习模型对多个样本进行m6A位点预测
#'
#' @param ml_fit 训练好的机器学习模型
#' @param feature_df 包含特征的数据框
#' @param positive_threshold 正样本阈值，默认为0.5
#' @return 包含预测概率和状态的数据框
#' @examples
#' \dontrun{
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' example_data <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#' prediction_multiple(rf_model, example_data, positive_threshold = 0.5)
#' }
#' @export
#' @import randomForest
#' @importFrom stats predict
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length",
                  "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in%
                  colnames(feature_df)))

  prediction_df <- feature_df

  prediction_df$RNA_type <- factor(prediction_df$RNA_type,
                                   levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  prediction_df$RNA_region <- factor(prediction_df$RNA_region,
                                     levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  dna_encoded <- dna_encoding(prediction_df$DNA_5mer)

  for(col in colnames(dna_encoded)) {
    dna_encoded[[col]] <- factor(dna_encoded[[col]], levels = c("A", "T", "C", "G"))
  }

  prediction_df_with_encoding <- cbind(prediction_df, dna_encoded)
  prediction_df_with_encoding$DNA_5mer <- NULL

  predicted_probs <- predict(ml_fit, newdata = prediction_df_with_encoding,
                             type = "prob")[, "Positive"]
  predicted_status <- ifelse(predicted_probs > positive_threshold,
                             "Positive", "Negative")

  feature_df$predicted_m6A_prob <- predicted_probs
  feature_df$predicted_m6A_status <- predicted_status

  return(feature_df)
}

#' 单样本m6A位点预测
#'
#' 该函数使用训练好的机器学习模型对单个样本进行m6A位点预测
#'
#' @param ml_fit 训练好的机器学习模型
#' @param gc_content GC含量数值
#' @param RNA_type RNA类型
#' @param RNA_region RNA区域
#' @param exon_length 外显子长度
#' @param distance_to_junction 到连接点的距离
#' @param evolutionary_conservation 进化保守性
#' @param DNA_5mer 5-mer DNA序列
#' @param positive_threshold 正样本阈值，默认为0.5
#' @return 包含预测概率和状态的命名向量
#' @examples
#' \dontrun{
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' prediction_single(rf_model, gc_content = 0.5, RNA_type = "mRNA",
#'                  RNA_region = "CDS", exon_length = 10, distance_to_junction = 8,
#'                  evolutionary_conservation = 0.5, DNA_5mer = "GGACA")
#' }
#' @export
#' @import randomForest
#' @importFrom stats predict
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5){
  single_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  single_df$RNA_type <- factor(single_df$RNA_type,
                               levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  single_df$RNA_region <- factor(single_df$RNA_region,
                                 levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  pred_result <- prediction_multiple(ml_fit, single_df, positive_threshold)

  returned_vector <- c(
    predicted_m6A_prob = pred_result$predicted_m6A_prob,
    predicted_m6A_status = pred_result$predicted_m6A_status
  )

  return(returned_vector)
}
