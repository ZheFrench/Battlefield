test_that("neighborhood of a spot given by identifier or coordinates", {
df <- expand.grid(x = 1:10, y = 1:10)
df$spot_id <- paste0("s", seq_len(nrow(df)))
df$cluster <- ifelse(df$x <= 5, "A", "B")
id <- df$spot_id[df$x == 5 & df$y == 5]
xy <- data.frame(x = 5, y = 5)

res <- get_neighborhood_spots(df, spot_id = id, k = 4)
expect_equal(nrow(res), 4)
expect_false(id %in% res$spot_id)
expect_true(all(abs(res$x - 5) + abs(res$y - 5) == 1))
expect_equal(get_neighborhood_spots(df, spot_id = xy, k = 4), res)

n <- count_neighborhood(df, spot_id = xy, k = 4)
expect_equal(n, count_neighborhood(df, spot_id = id, k = 4))
expect_equal(unique(n$cluster), id)
expect_equal(n$count[match(c("A", "B"), n$neighborhood)], c(3, 1))

expect_error(get_neighborhood_spots(df, spot_id = data.frame(x = 5.4, y = 5)),
    "No spot found")
expect_error(get_neighborhood_spots(df, spot_id = "nope"), "No spot found")
expect_error(get_neighborhood_spots(df, spot_id = c("s1", "s2")), "single")
})
