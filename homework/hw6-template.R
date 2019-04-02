fn = "mutt_gamete_event.tsv"
data = read.csv(fn, header=TRUE, sep="\t")
events = data$event
distance = data$num_base_pairs_to_prev

lnlf = function(r) {
  lnl = 0.0;
  for (i in 1:length(events)) {
    lnl.for.event = 0;
    if (as.character(events[i]) == 'E') {
      et = "end of chromosome";
    } else {
      et = "Recombination";
    }
    d = distance[i];
    print(paste("Calc lnl.for.event here for ", et, "at", as.character(d), "bases after previous event |  r=", r))
    lnl = lnl + lnl.for.event;
  }
  return(lnl);
}
lnlf(0.1)
num.opt.answer = optimize(lnlf, interval = c(0, 1), maximum = TRUE)

mle.r = num.opt.answer$maximum
crit.lnL = num.opt.answer$objective - 1.92
