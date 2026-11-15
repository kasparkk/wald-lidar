# Called by treeflow: segment trees with lidR and hand back point -> tree ids.
# Args: input.las  out_dir  algorithm(dalponte|silva|li2012)  ws_a  ws_b  hmin  chm_res
suppressPackageStartupMessages(library(lidR))
a <- commandArgs(TRUE)
inp <- a[1]; out <- a[2]; algo <- a[3]
ws_a <- as.numeric(a[4]); ws_b <- as.numeric(a[5]); hmin <- as.numeric(a[6]); res <- as.numeric(a[7])
set_lidr_threads(2)
las <- readLAS(inp, select = "xyzcr0")
if (sum(las$Classification == 2L) < 0.01 * npoints(las)) las <- classify_ground(las, csf())
nlas <- normalize_height(las, tin())
chm <- rasterize_canopy(nlas, res, pitfree(subcircle = 0.15))
ws <- function(x) pmax(2, ws_a + ws_b * x)            # window diameter in metres, grows with height
ttops <- locate_trees(chm, lmf(ws, hmin = hmin))
seg <- switch(algo,
  dalponte = dalponte2016(chm, ttops),
  silva    = silva2016(chm, ttops),
  li2012   = li2012())
nlas <- segment_trees(nlas, seg)
pid <- nlas@data[["pid"]]
if (is.null(pid)) stop("point id attribute 'pid' missing")
tid <- nlas$treeID; tid[is.na(tid)] <- 0L
data.table::fwrite(data.table::data.table(pid = pid, tree = tid), file.path(out, "lidr_points.csv"))
if (algo != "li2012") {
  xy <- sf::st_coordinates(ttops)
  data.table::fwrite(data.table::data.table(tree = ttops$treeID, x = xy[, 1], y = xy[, 2], z = xy[, 3]),
                     file.path(out, "lidr_tops.csv"))
}
cat("lidR done:", length(unique(tid[tid > 0])), "trees\n")
