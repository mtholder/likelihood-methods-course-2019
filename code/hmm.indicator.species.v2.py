# v2 : estimate habitat transition rate

from math import exp,log
from scipy import optimize

def calc_t_rate(r_rates): 
	def scipy_ln_like0(x):
		llx,fprbs,rprbs=foward_backward(observations,states,start_probability,x)
		return -llx

	bounds = [ (0.001,0.999) ]
	best, val, d = optimize.fmin_l_bfgs_b(scipy_ln_like0, r_rates, approx_grad=True, bounds=bounds)
	solution = list(best)
	ln_l = -scipy_ln_like0(solution)
	solution.append(ln_l)
	
	return solution


def emission_probability(truth,observed):
	if truth=='habitat A':
		if observed=='00':
			return 0.1
		elif observed=='10':
			return 0.7
		elif observed=='01':
			return 0.01
		elif observed=='11':
			return 0.19
	elif truth=='habitat B':
		if observed=='00':
			return 0.2
		elif observed=='10':
			return 0.02
		elif observed=='01':
			return 0.63
		elif observed=='11':
			return 0.15
	return

def foward_backward(obs, states, start_p, tprob):


	alpha=[{} for j in range(len(obs))] # forward:: alpha[j][X] is probability that true state is X at sample j (starts at 0)
	beta= [{} for j in range(len(obs))] # backward:: beta[j][X] is probability that true state is X at sample j (starts at 0)

# moved this into function
	transition_probability=[{} for j in range(Samples-1)] # 
	for x1 in range(Samples-1): # 
		transition_probability[x1] ={ 'habitat A': {'habitat A':1-tprob[0],'habitat B':tprob[0]}, 'habitat B': {'habitat A':tprob[0],'habitat B':1-tprob[0]} }



	lnFactor=0.0

	for y in states:
		alpha[0][y] = start_p[y] * emission_probability(y,obs[0])
		beta[len(obs)-1][y] = 1.0


	for t in range(1, len(obs)):
		for y in states:
			alpha[t][y] = 0.0
			for y0 in states: # y0 is state at t-1
				alpha[t][y] +=alpha[t-1][y0] * transition_probability[t-1][y0][y] * emission_probability(y,obs[t])

		normalizer = max(alpha[t]['habitat A'],alpha[t]['habitat B'])
		lnFactor+=log(normalizer)

		for y in states:
			alpha[t][y] = alpha[t][y]/normalizer

	LLobs=lnFactor+log(alpha[len(obs)-1]['habitat A']+alpha[len(obs)-1]['habitat B'])	

	for t in range(len(obs)-2,-1,-1):
		for y in states:
			beta[t][y] = 0.0 # y is state at t
			for y0 in states: # y0 is state at t+1
				beta[t][y] +=beta[t+1][y0] * transition_probability[t][y][y0] * emission_probability(y0,obs[t+1])

	return LLobs,alpha,beta




# main program

observations=['01','01','00','01','01','01','01','01','00','11','10','11','11','10','11','11','10','11','01','00','10','10','10','00','10']
Samples=len(observations)

states = ('habitat A','habitat B')
start_probability = {'habitat A':0.5,'habitat B':0.5}






t_rate=[0.5] # habitat transition rate to be estimated
zsol= calc_t_rate(t_rate) 
tr = zsol[0]

llx,fprbs,rprbs=foward_backward(observations,states,start_probability,[tr])

print("LL = ",llx," mle Transition rate ",tr)

postProb=[{} for j in range(len(fprbs))] # 
for j in range(len(fprbs)):
	denom=0.0
	for y in states: 
		denom+=(fprbs[j][y]*rprbs[j][y])
	for y in states: 
		postProb[j][y]=(fprbs[j][y]*rprbs[j][y])/denom


for j in range(Samples):
	print(j,observations[j],postProb[j]['habitat A'],postProb[j]['habitat B'])






