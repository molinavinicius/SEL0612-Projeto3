from inspect import isfunction
from math import sqrt, sin, pi

mi0 = 4e-7 * pi
eps0 = 8.854187817e-12

def fdtd_gen(i_max, delta, delta_t, coefs, f_meio=None, j_max=None,):
	""" Cria um gerador para a simulação FDTD
	Este gerador retorna a cada passo uma iteração no tempo (n) da simulação FDTD
	"""
	# parametro opcional. se não fornecido, malha quadrada
	if j_max == None:
		j_max = i_max 
	
	#Inicializar as variaveis
	Ez_ant = Hy_ant = Hx_ant = Ez = Hx = Hy = Jz_fonte = Mx_fonte = My_fonte = [[0 for i in range(0, i_max)] for j in range(0, j_max)]
	
	n = 0
	#Loop de iteração
	while True:
		yield Ez, Hx, Hy

		n = n + 1

		# avançar um passo
		Ez_ant = Ez
		Hx_ant = Hx
		Hy_ant = Hy

		#calcula Ez
		Ez = (
			[[0 for j in range (0, j_max)]]
			+ [[0] + [passo_Ez(i, j, Ez_ant, Hy_ant, Hx_ant, Jz_fonte, delta, coefs) for j in range(1, j_max - 1)] + [0] for i in range(1, i_max - 1)]
			+ [[0 for j in range(0, j_max)]]
		)

		#se fornecida uma excitação, aplicá-la no ponto central
		if isfunction(f_meio):
			Ez[round(i_max/2)][round(j_max/2)] = f_meio(n*delta_t)

		#calcular Hx e Hy
		Hx = [[0] + [passo_Hx(i, j, Hx_ant, Ez, Mx_fonte, delta, coefs) for j in range(1, j_max)] for i in range(0, i_max)]

		Hy = (
			[[passo_Hy(i, j, Hy_ant, Ez, My_fonte, delta, coefs) for j in range(0, j_max)] for i in range(0, i_max - 1)]
			+ [[0 for j in range(0, j_max)]]
		)

def passo_Ez(i, j, Ez_ant, Hy_ant, Hx_ant, Jz_fonte, delta, coefs):
	return coefs[i][j]['Ca'] * Ez_ant[i][j] + coefs[i][j]['Cb'] * (Hy_ant[i][j] - Hy_ant[i-1][j] + Hx_ant[i][j] - Hx_ant[i][j+1] - Jz_fonte[i][j]*delta)

def passo_Hx(i, j, Hx_ant, Ez, Mx_fonte, delta, coefs):
	return coefs[i][j]['Da']*Hx_ant[i][j] + coefs[i][j]['Db'] * (Ez[i][j-1] - Ez[i][j] - Mx_fonte[i][j]*delta)

def passo_Hy(i, j, Hy_ant, Ez, My_fonte, delta, coefs):
	return coefs[i][j]['Da']*Hy_ant[i][j] + coefs[i][j]['Db'] * (Ez[i+1][j] - Ez[i][j] - My_fonte[i][j]*delta)
