#from distutils.core import setup # Finish this
setup(
  name = 'snp_haplotype',         
  packages = ['snp_haplotype'],  
  version = '0.1',      # Start with a small number and increase it with every change you make
  license='MIT',        
  description = 'Finds SNP haplotypes between samples and produces a Minimum Spanning Tree showing their relationship',
  author = 'Zac Scurlock',                   
  author_email = 'zac.scur@hotmail.com',      
  url = 'https://github.com/ZScurlock/SNP_haplotype',  
  download_url = 'https://github.com/user/reponame/archive/v_01.tar.gz',    # I explain this later on
  keywords = ['SNP', 'Haplotype', 'Haplotype network', 'Minimum spanning tree'],
  install_requires=[            
          'pandas',
          'numpy',
          'networkx'
          'matplotlib'
          'pairsnp'
    
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Bioinformaticians',      
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.6.9',
  ],
)
