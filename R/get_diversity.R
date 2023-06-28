get_diversity <- function(data_source, round = TRUE, sel_method = "") {
  # helper function
  get_diversity_taxonomic <- function(data_matrix, sample_size) {
    hill0 <- function(data, sample_size) {
      data_sub <- data[data > 0]
      data_sum <- sum(data_sub)
      if (sample_size <= data_sum) {
        res <- sum(1 - exp(lchoose(data_sum - data_sub, sample_size) -
          lchoose(data_sum, sample_size)))
        return(res)
      } else {
        return(0)
      }
    }
    fk_hat <- function(data, sample_size) {
      data_sub <- data[data > 0]
      data_sum <- sum(data_sub)
      if (sample_size <= data_sum) {
        sub <- function(k) {
          sum(exp(lchoose(data_sub, k) + lchoose(data_sum -
            data_sub, sample_size - k) - lchoose(
            data_sum,
            sample_size
          )))
        }
        res <- sapply(1:sample_size, sub)
        return(res)
      } else {
        return(0)
      }
    }
    hill1 <- function(data, sample_size) {
      data_sub <- data[data > 0]
      data_sum <- sum(data_sub)
      if (sample_size <= data_sum) {
        k <- 1:sample_size
        res <- exp(-sum(k / sample_size * log(k / sample_size) *
          fk_hat(data_sub, sample_size)))
        return(res)
      } else {
        return(0)
      }
    }
    hill2 <- function(data, sample_size) {
      data_sub <- data[data > 0]
      data_sum <- sum(data_sub)
      if (sample_size <= data_sum) {
        res <- 1 / (1 / sample_size + (1 - 1 / sample_size) * sum(data_sub *
          (data_sub - 1) / data_sum / (data_sum - 1)))
        return(res)
      } else {
        return(0)
      }
    }
    est_n_0 <-
      sapply(sample_size, function(n) {
        apply(data_matrix, 1, hill0, sample_size = n)
      })
    est_n_1 <-
      sapply(sample_size, function(n) {
        apply(data_matrix, 1, hill1, sample_size = n)
      })
    est_n_2 <-
      sapply(sample_size, function(n) {
        apply(data_matrix, 1, hill2, sample_size = n)
      })
    hill_diversity <-
      cbind(est_n_0, est_n_1, est_n_2, est_n_1 -
        est_n_2, est_n_2 / est_n_1, est_n_1 / est_n_0) %>%
      as.data.frame() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(sample_id = row.names(data_matrix)) %>%
      tibble::column_to_rownames("sample_id") %>%
      purrr::set_names(
        nm = c(
          "n0",
          "n1", "n2", "n1_minus_n2", "n2_divided_by_n1", "n1_divided_by_n0"
        )
      )
    return(hill_diversity)
  }

  if (
    isTRUE(round)
  ) {
    data_matrix <-
      data_source %>%
      tibble::column_to_rownames("sample_id") %>%
      dplyr::mutate_all(., .f = floor) %>%
      as.matrix() %>%
      round()
  } else {
    data_matrix <-
      data_source %>%
      tibble::column_to_rownames("sample_id") %>%
      as.matrix()
  }

  sample_size <-
    apply(data_matrix, 1, sum) %>%
    floor() %>%
    min()

  div <-
    get_diversity_taxonomic(
      data_matrix = data_matrix,
      sample_size = sample_size
    )

  res <-
    div %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::relocate(sample_id) %>%
    as.data.frame() %>%
    tibble::as_tibble()
  return(res)
}
