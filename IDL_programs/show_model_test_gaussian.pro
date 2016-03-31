@read_params
; show results obtained using the model_Test_Gaussian
pro show_model_test_gaussian, file_data_obs, dir_outputs

dir_outputs='~/Dropbox/Temporary/Cpp-Playground/Data_Tests/outputs/'
file_data_obs='~/Dropbox/Temporary/Cpp-Playground/Data_Tests/test-gauss.txt'

data_struc=read_obsdata(file_data_obs)

file_prop_covarmat=dir_outputs + 'proposals_covarmats_chain-0.txt'
file_prop_moves=dir_outputs + 'proposals_moves.txt'
file_prop_sigma=dir_outputs + 'proposals_sigmas.txt'
file_prop_mu=dir_outputs + 'proposals_mus_chain-0.txt'

file_stats=dir_outputs + 'stat_criteria.txt'
file_parallel_temp=dir_outputs + 'parallel_tempering.txt'
file_params=dir_outputs + 'params_chain-0.txt'

;print, 'Reading the covariance matrices...'
;covarmats_struc=read_proposal_covarmat(file_prop_covarmat)
print, 'Reading Pmove and moved variables...'
moves_struc=read_proposal_moves(file_prop_moves)
print, 'Reading the sigma scaling factors...'
sigma_struc=read_proposal_sigmas(file_prop_sigma)
print, 'Reading the averaged target variables (mu)...'
mu_struc=read_proposal_mus(file_prop_mu)
print, 'Reading the logLikelihood, logPrior, logPosterior...'
stats_struc=read_stat_criteria(file_stats)
;print, 'Reading the parallel tempering parameters...'
;parallel_temp_struc=read_parallel_tempering(file_parallel_temp)
print, 'Reading the variables and constants of the model...'
params_struc=read_params(file_params)

Nclass=20
Nsubdata=moves_struc.Ndata/Nclass;
acceptance=dblarr(moves_struc.Nchains, Nclass);
xiter=dblarr(moves_struc.Nchains, Nclass)
for i=0, Nclass-1 do begin
	for chain=0, moves_struc.Nchains-1 do begin
		moved=where(moves_struc.moved[chain,i*Nsubdata:(i+1)*Nsubdata-1] eq 1)
		notmoved=where(moves_struc.moved[chain,i*Nsubdata:(i+1)*Nsubdata -1] eq 0)
		acceptance[chain,i]=1. * n_elements(moved)/Nsubdata
	endfor
	xiter[i]=i*Nsubdata + 1.*Nsubdata/2
endfor
plot, xiter, 1 - acceptance[0,*], /xst

if params_struc.cons_names[0] eq 'None' then begin
	params=params_struc.vars; Only true if you have Ncons=0... otherwise, you need to reconstruct the vector
endif else begin
	print, 'Warning: The function read_params() of read_params.pro does not yet read relax!'
	print, 'Need to be implemented and tested'
	stop
	pcons=where(params_struc.relax eq 0)
	pvars=where(params_struc.relax eq 1)
	params=dblarr(total(pcons) + total(pvars), params_struc.Ndata)
	params[pcons,*]=params_struc.cons
	params[pvars,*]=params_struc.vars
endelse

meds=dblarr(params_struc.Nvars)
for i=long(0), params_struc.Nvars-1 do begin
	meds[i]=median(params[i,*])
endfor

model=model_test_gaussian(data_struc.vars[0,*], meds)		

plot, data_struc.vars[0,*], data_struc.vars[1,*]
oplot, data_struc.vars[0,*], model, color=fsc_color('blue')

stop
end

; Identical model to the one used by models.cpp and names 'model_Test_Gaussian'
function model_test_gaussian, x, params

	model=params[0] * exp( -0.5 * ( x - params[2])^2 / params[1]^2) + params[3]
	
return, model
end