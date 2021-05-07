using Revise
using MCMCChains
c1 = Chains(randn(100,5, 1), info = (start_time=time(), stop_time = time()+1))
c2 = Chains(randn(100,5, 1), info = (start_time=time(), stop_time = time()+1))
c = chainscat(c1, c2)

c1 = Chains(randn(100,5, 1), info = (start_time=[time()], stop_time = [time()+1]))
c2 = Chains(randn(100,5, 1), info = (start_time=[time()], stop_time = [time()+1]))
c = chainscat(c1, c2)


MCMCChains.compute_duration(c)