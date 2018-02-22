VirusSeqPipeline: A VirusSeq pipeline for processing multiple RNAseq samples: a command line approach 
=========

The VirusSeq software is a next generation sequencing software of human cancer tissue. Its goal is to identify viruses and their integration sites in human cancer cells. VirusSeq has originally been developped by [Chen and collaborators](https://academic.oup.com/bioinformatics/article/29/2/266/202055 "VirusSeq: software to identify viruses and their integration sites using next-generation sequencing of human cancer tissue") and reported in a 2012 paper published in Bioinformatics. Although the software is open source and freely available, its installation is not at the ease of the novice user. VirusSeq has not been designed to effectively and efficiently process multiple samples. In its present shape, it requires the user to write a shell script to every single sample to process. Furthermore, the user has to engage in the tedious and error prone tasks of path modification. All these limitations of the VirusSeq software may tremendously limit its use in the scientific community. 

In this project, we aim at addressing the aforementioned limitations by designing an easy to install, easy to implement workflow for the VirusSeq software. Our approach does not require the user to have knowledge of any shell command. In this user guide documentation, we present the file organization and a more compact package for the user to install on his computer. 



VirusSeqPipeline description 
=========
This pipeline is composed of the MOSAIK aligner, reference genome and annotation files and PERL programs. The PERL programs are originally from the VirusSeq software. In the full version of this workflow, we provided the user with the reference genome and annotation files made available by the authors of VirusSeq and available at : http://odin.mdacc.tmc.edu/~xsu1/VirusSeq.html. Also, see the original user guide at the same URL for more information. 

The file organization of our workflow is presented as follow: 

- the directory VirusSeq/Mosaik_bin contains all the Mosaik aligner executable tools 

- the directory VirusSeq/MosaikJumpDb contains the following files: 

- the reference genome and annotation files: 

1. hg19.fa  

2. gibVirus.fa  

3. hg19Virus.fa  

4. hg19Virus_refGene_RIS.txt  

5. Spanner_anchor_hg19Virus.txt 

- a SHELL script Jump_file_builder.sh to generate the jump files from the reference genome and annotation files 

- the VirusSeq/VirusSeq_Script directory contains the following PERL programs: 

1. Spanner_cross_converter.pl 

2. VirusSeq_Detection.pl 

3. VirusSeq_Integration.pl 

- the VirusSeq/samples directory is where the user is supposed to place the samples to be processed. In section 2, we provide detailed direction on how the user is expected to organize his FASTQ samples. 

- the VirusSeq/results directory contains the files resulting from the processing of each sample organized in folders named after each sample. 

- The VirusSeq folder contains the file runAll.sh which is a SHELL script containing the command to run the VirusSeq pipeline on all the samples contained in the VirusSeq/samples directory. This directory also contain other archive files for the MOSAIK aligner, the reference genome and annotations file as well as a copy of the original VirusSeq package. In addition, we also provide a copy of the eagleView software used to preview the viruses sequences in .ACE format. On some computers and for users with limited privileges, the Mosaik tools might not be allow to run as executable programs. This usually happens when users attempt to run this workflow from an external media. Knowing that this may limit the use of this workflow, we provide the file AllowExec.txt to help users figure out how to allow the Mosaik toole to run as executable. 


We bring to the user's attention to the fact that altering or changing this file organization may prevent VirusSeq from running correctly. 



VirusSeqPipeline data input
=========
This pipeline requires the user to input RNA-seq data in a specific compressed FASTQ format (.fq.gz) with paired-end reads. If for example the user is to process a sample "L526401A", the FASTQ files "L526401A_1.fq.qz" and "L526401A_2.fq.gz" are required. We expect the user to create a folder named "L526401A" containing the files "L526401A_1.fq.qz" and "L526401A_2.fq.gz". The folder "L526401A" is later to be placed in the directory VirusSeq/samples for its processing. 

In general, Each folder in the VirusSeq/samples directory is a sample and a sample is made of two FASTQ files with the respective extensions "_1.fq.gz" and "_2.fq.gz". The folder in the VirusSeq/samples directory must be named after the name of the sample it represents. 

It is important to note that after the processing of each sample, two additional folders named respectively "Gig" and "SV_gDNA" are created inside in subfolder of the VirusSeq/samples directory. The "Gig" folder contains the sequence files (in .ACE format) of all the viruses reference found for that sample. Each .ACE file can be preview using the freely availabel EagleView package we provided along with this workflow. EagleView is also available for download at: http://www.niehs.nih.gov/research/resources/software/biostatistics/eagleview/index.cfm. 

The "SV_gDNA" folder received the results from the virus integration sites. 



How to install and run VirusSeqPipeline
=========
Please, follow the instructions below to set up the workflow and run it: 

- Clone VirusSeqPipeline

	git clone https://github.com/rosericazondekon/VirusSeqPipeline.git

- Now open a shell terminal (CTRL+TAB+T on Ubuntu) and change your directory (using the "CD" command) to your VirusSeqPipeline directory 

- Make sure that all the Mosaik aligner tools are allowed to be executed as files. To do this, navigate to the VirusSeq/Mosaik_bin directory and right click on each file, then select "Properties", next select "Permissions" and tick the box "Allow executing file as program". If for some reasons, this is impossible, follow the instructions in "AllowExec.txt" file. 

- Put your samples organized as explained in section 2 in the VirusSeq/samples directory. Do remember that each sample directory must contain two compressed FASTQ files ending respectively by the extensions "_1.fq.gz" and "_2.fq.gz" 

- And run the following shell command:
```shell
bash runAll.sh
```

Do notice that depending on the number of samples, you might want to make sure that hundreds of Gigabytes of space are available on your media. The entire processing may therefore take several hours or days to complete.






You simply download/clone the source files into the directory you keep your python scripts. 
Then you add that directory to your PYTHONPATH.
In Unix/Linux you can do this by adding into your .bashrc file the two following lines:
	
	export PYTHONPATH=/home/your_username/your_python_libs:$PYTHONPATH
	export PYTHONPATH=/home/your_username/your_python_libs/multinetx:$PYTHONPATH


After request of some users I give here an example of how to "install" multinetx:

Create a directory for your python libraries (if you do not have already)

	mkdir your_python_libs

Enter this directory

	cd your_python_libs

Clone the multinetx

	git clone https://github.com/nkoub/multinetx.git

Add multinetx to your PYTHONPATH by adding to your .bashrc the two following lines:

	export PYTHONPATH=/home/your_username/your_python_libs:$PYTHONPATH
	export PYTHONPATH=/home/your_username/your_python_libs/multinetx:$PYTHONPATH



How to use multiNetX
=========


#### Import standard libraries for numerics


```python
import numpy as np
```

#### Import the package MultiNetX

```python
import multinetx as mx
```

# Create a multiplex 1st way

#### Create three Erd"os- R'enyi networks with N nodes for each layer

```python
N = 5
g1 = mx.generators.erdos_renyi_graph(N,0.5,seed=218)
g2 = mx.generators.erdos_renyi_graph(N,0.6,seed=211)
g3 = mx.generators.erdos_renyi_graph(N,0.7,seed=208)
```
#### Create an 3Nx3N lil sparse matrix. It will be used to describe the layers interconnection

```python
adj_block = mx.lil_matrix(np.zeros((N*3,N*3)))
```
#### Define the type of interconnection among the layers (here we use identity matrices thus connecting one-to-one the nodes among layers)

```python
adj_block[0:  N,  N:2*N] = np.identity(N)    # L_12
adj_block[0:  N,2*N:3*N] = np.identity(N)    # L_13
adj_block[N:2*N,2*N:3*N] = np.identity(N)    # L_23
    
# use symmetric inter-adjacency matrix
adj_block += adj_block.T
```
#### Create an instance of the MultilayerGraph class

```python
mg = mx.MultilayerGraph(list_of_layers=[g1,g2,g3],
                        inter_adjacency_matrix=adj_block)
```
#### Weights can be added to the edges

```python
mg.set_edges_weights(intra_layer_edges_weight=2,
                     inter_layer_edges_weight=3)
```
# Create a multiplex 2nd way

```python
mg = mx.MultilayerGraph()
```
#### Add layers

```python
mg.add_layer(mx.generators.erdos_renyi_graph(N,0.5,seed=218))
mg.add_layer(mx.generators.erdos_renyi_graph(N,0.6,seed=211))
mg.add_layer(mx.generators.erdos_renyi_graph(N,0.7,seed=208))
```
#### Create an instance of the MultilayerGraph class

```python
mg.layers_interconnect(inter_adjacency_matrix=adj_block)
```
#### Weights can be added to the edges

```python
mg.set_edges_weights(intra_layer_edges_weight=2,
                     inter_layer_edges_weight=3)
```

The object mg inherits all properties from Graph of networkX, so that
we can calculate adjacency or Laplacian matrices, their eigenvalues, etc.




How to plot multiplex networks
=========
     

##### Import standard libraries

```python
import numpy as np
import matplotlib.pyplot as plt
```
##### Import the package MultiNetX

```python
import multinetx as mx
```
##### Create three Erd"os- R'enyi networks with N nodes for each layer

```python
N = 50
g1 = mx.erdos_renyi_graph(N,0.07,seed=218)
g2 = mx.erdos_renyi_graph(N,0.07,seed=211)
g3 = mx.erdos_renyi_graph(N,0.07,seed=208)
```
### Edge colored nertwork (no inter-connected layers)

##### Create the multiplex network

```python
mg = mx.MultilayerGraph(list_of_layers=[g1,g2,g3])
```
##### Set weights to the edges

```python
mg.set_intra_edges_weights(layer=0,weight=1)
mg.set_intra_edges_weights(layer=1,weight=2)
mg.set_intra_edges_weights(layer=2,weight=3)
```
##### Plot the adjacency matrix and the multiplex networks

```python
fig = plt.figure(figsize=(15,5))
ax1 = fig.add_subplot(121)
ax1.imshow(mx.adjacency_matrix(mg,weight='weight').todense(),
		  origin='upper',interpolation='nearest',cmap=plt.cm.jet_r)
ax1.set_title('supra adjacency matrix')

ax2 = fig.add_subplot(122)
ax2.axis('off')
ax2.set_title('edge colored network')
pos = mx.get_position(mg,mx.fruchterman_reingold_layout(g1),
					  layer_vertical_shift=0.2,
					  layer_horizontal_shift=0.0,
					  proj_angle=47)
mx.draw_networkx(mg,pos=pos,ax=ax2,node_size=50,with_labels=False,
				 edge_color=[mg[a][b]['weight'] for a,b in mg.edges()],
				 edge_cmap=plt.cm.jet_r)
plt.show()
```

![png](plot_multiplex_networks_files/plot_multiplex_networks_12_0.png)


### Regular interconnected multiplex

##### Define the type of interconnection between the layers

```python
adj_block = mx.lil_matrix(np.zeros((N*3,N*3)))

adj_block[0:  N,  N:2*N] = np.identity(N)    # L_12
adj_block[0:  N,2*N:3*N] = np.identity(N)    # L_13
#adj_block[N:2*N,2*N:3*N] = np.identity(N)    # L_23
adj_block += adj_block.T
```
##### Create an instance of the MultilayerGraph class

```python
mg = mx.MultilayerGraph(list_of_layers=[g1,g2,g3], 
						inter_adjacency_matrix=adj_block)

mg.set_edges_weights(inter_layer_edges_weight=4)

mg.set_intra_edges_weights(layer=0,weight=1)
mg.set_intra_edges_weights(layer=1,weight=2)
mg.set_intra_edges_weights(layer=2,weight=3)
```
##### Plot the adjacency matrix and the multiplex networks

```python
fig = plt.figure(figsize=(15,5))
ax1 = fig.add_subplot(121)
ax1.imshow(mx.adjacency_matrix(mg,weight='weight').todense(),
		  origin='upper',interpolation='nearest',cmap=plt.cm.jet_r)
ax1.set_title('supra adjacency matrix')

ax2 = fig.add_subplot(122)
ax2.axis('off')
ax2.set_title('regular interconnected network')
pos = mx.get_position(mg,mx.fruchterman_reingold_layout(mg.get_layer(0)),
					  layer_vertical_shift=1.4,
					  layer_horizontal_shift=0.0,
					  proj_angle=7)
mx.draw_networkx(mg,pos=pos,ax=ax2,node_size=50,with_labels=False,
				 edge_color=[mg[a][b]['weight'] for a,b in mg.edges()],
				 edge_cmap=plt.cm.jet_r)
plt.show()
```

![png](plot_multiplex_networks_files/plot_multiplex_networks_19_0.png)


### General multiplex multiplex 

##### Define the type of interconnection between the layers

```python
adj_block = mx.lil_matrix(np.zeros((N*4,N*4)))

adj_block[0  :  N ,   N:2*N] = np.identity(N)   # L_12
adj_block[0  :  N , 2*N:3*N] = np.random.poisson(0.005,size=(N,N))   # L_13
adj_block[0  :  N , 3*N:4*N] = np.random.poisson(0.006,size=(N,N))   # L_34
adj_block[3*N:4*N , 2*N:3*N] = np.random.poisson(0.008,size=(N,N))   # L_14
adj_block += adj_block.T
adj_block[adj_block>1] = 1
```
##### Create an instance of the MultilayerGraph class

```python
mg = mx.MultilayerGraph(list_of_layers=[g1,g2,g3,g1],
						inter_adjacency_matrix=adj_block)

mg.set_edges_weights(inter_layer_edges_weight=5)

mg.set_intra_edges_weights(layer=0,weight=1)
mg.set_intra_edges_weights(layer=1,weight=2)
mg.set_intra_edges_weights(layer=2,weight=3)
mg.set_intra_edges_weights(layer=3,weight=4)
```
##### Plot the adjacency matrix and the multiplex networks

```python
fig = plt.figure(figsize=(15,5))
ax1 = fig.add_subplot(121)
ax1.imshow(mx.adjacency_matrix(mg,weight='weight').todense(),
		  origin='upper',interpolation='nearest',cmap=plt.cm.jet_r)
ax1.set_title('supra adjacency matrix')

ax2 = fig.add_subplot(122)
ax2.axis('off')
ax2.set_title('general multiplex network')
pos = mx.get_position(mg,mx.fruchterman_reingold_layout(mg.get_layer(0)),
					  layer_vertical_shift=.3,
					  layer_horizontal_shift=0.9,
					  proj_angle=.2)
mx.draw_networkx(mg,pos=pos,ax=ax2,node_size=50,with_labels=False,
				 edge_color=[mg[a][b]['weight'] for a,b in mg.edges()],
				 edge_cmap=plt.cm.jet_r)
plt.show()
```

![png](plot_multiplex_networks_files/plot_multiplex_networks_26_0.png)



    




Copyright
=========

(C) Copyright 2013-2015, Nikos E Kouvaris

Each file in this folder is part of the multiNetX package.

multiNetX is part of the deliverables of the LASAGNE project 
(multi-LAyer SpAtiotemporal Generalized NEtworks),
EU/FP7-2012-STREP-318132 (http://complex.ffn.ub.es/~lasagne/)

multiNetX is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

multiNetX is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.