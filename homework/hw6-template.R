fn = "mutt_gamete_event.tsv"
data = read.csv(fn, header=TRUE, sep="\t")
events = data$event
distance = data$num_base_pairs_to_prev
num.chrom = sum(as.character(events) == "E")
num.crossovers = sum(as.character(events) == "R")

lnlf = function(r) {
  lnl = 0.0;
  prev.was.recomb = FALSE;
  for (i in 1:length(events)) {
    d = distance[i];
    lnl.for.event = 0;
    if (prev.was.recomb) {
        print("Add something to lnl.for.event following recomb");
    } else {
        print("Add something to lnl.for.event first event on chromosome");
    }
    if (as.character(events[i]) == 'E') {
      print("Perhaps add something to lnl.for.event for no recomb at this point");
      prev.was.recomb = FALSE;
    } else {
      print("Perhaps add something to lnl.for.event for recomb at this point");
      prev.was.recomb = TRUE;
    }
    lnl = lnl + lnl.for.event;
  }
  return(lnl);
}
lnlf(0.1)
num.opt.answer = optimize(lnlf, interval = c(0, 1), maximum = TRUE)

mle.r = num.opt.answer$maximum
crit.lnL = num.opt.answer$objective - 1.92
