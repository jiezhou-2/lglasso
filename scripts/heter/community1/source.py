import os,time,ConfigParser
def main():
		
		for m  in ['20','40']:
			for n in ['30', '10','20']:
				for e in ['0','0.5','1']:
					for p in ['80']:
						identifier='subject='+m+'_time='+n+'_coef='+e+'_nodes='+p
						outfiles =  identifier + '.Rd'
						RFileName = 'heter'+identifier  +'.R' 
						RFile = open(RFileName, 'w')
						RFile.write('library("lglasso")\n')
						RFile.write('library("glasso")\n')
						RFile.write('library("BDgraph")\n')
						RFile.write('library("GGMselect")\n')
						RFile.write('m='+m+'\n')
						RFile.write('n='+n+'\n')
						RFile.write('e='+e+'\n')
						RFile.write('p='+p+'\n')					
						RFile.write('Nsim=50\n')
						RFile.write('rho=vector("list",5)\n') 
						RFile.write('rho[[1]]=seq(0.001,0.3,length=20)\n')
						RFile.write('rho[[2]]=seq(0.001,0.3,length=20)\n')
						RFile.write('rho[[3]]=seq(0.001,0.3,length=20)\n')
						RFile.write('rho[[4]]=seq(0.00001,0.1,length=10)\n')
						RFile.write('rho[[5]]=seq(0.00001,0.1,length=10)\n')
						RFile.write('set.seed(m'+ '+n'+'+e'+'+p'+')'+'\n')                 			
						RFile.write('simures=power_compare1(m=m,n=n,p=p,coe=c(2,e,e),l=Nsim,rho=rho,prob=0.01,heter=T)\n')
						RFile.write('save(simures,file="'+outfiles +'")\n')
						RFile.write('q("no")\n')
						RFile.close()
						shFileName = 'heter'+identifier  +'.sh'
						shFile = open(shFileName, 'w')
						shFile.write('#!/bin/bash\n')
						shFile.write('#SBATCH --job-name=lglasso'+m+'_'+n+'_'+e+'_'+p+'\n')
						shFile.write('#SBATCH --nodes=1\n')
						shFile.write('#SBATCH --ntasks-per-node=1\n')
						shFile.write('#SBATCH --time=72:00:00\n')
						shFile.write('#SBATCH --mail-type=BEGIN,END,FAIL\n')
						shFile.write('#SBATCH --output=real.out\n')
						shFile.write('#SBATCH --error=real.err\n')
						shFile.write('time R CMD BATCH '  + RFileName + ' a' +  identifier  +'.out\n')
						shFile.close()
						os.system ('sbatch ' + shFileName)
						#time.sleep(1) #delay 1/10 of a second between job submissions

if __name__=="__main__":
	main()
	
