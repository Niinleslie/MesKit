# Schism2Fishplot

## Data Requirement

1. Cellularity of each element in the fraction table of fishplot should be larger than 0. 
2. The total cellularities of clusters in the same nest should be smaller than 100. 
3. According to Lineage Precedence Rule, the sum of cellularities from child clusters should not be larger than their parent cluster. This depends on the evolution relationship referred by SCHISM and the statue of your data.  
4. To complement upper "Prone-to-error" rules, you should obey other restrictions in cancer  sub-clone evolution research which could help your data to process successfully.



## Usage

1. confirm the directory of the cluster.tsv and loci.tsv

   ```
   dir.cluster.tsv
   dir.loci.tsv
   ```

2. use the `prepareSchismInput` function in the Schism2Fishplot.R and set the output directory

   ```R
   prepareSchismInput(dir.cluster.tsv, dir.loci.tsv, dir.output)
   ```

   input: 

   - `hu.cluster.tsv`
   - `hu.loci.tsv`

   output: 

   - `W.clusterEstimates.tsv`
   - ` W.mutation-to-cluster.tsv`

3. check the W.yaml in Schism folder and make sure the`working_dir`、`mutation_to_cluster_assignment`、`mutation_cellularity_input`、`output_prefix` for running Schism is right.

   ```yaml
   
   working_dir: /home/ninomoriaty/R_Project/EvolCancer/EvolCancer/Schism
   
   mutation_to_cluster_assignment: W.mutation-to-cluster.tsv
   
   mutation_cellularity_input: W.clusterEstimates.tsv
   
   output_prefix: 0_W_results
   
   cellularity_estimation: other
   
   hypothesis_test:
     test_level: clusters
     significance_level: 0.05
     store_pvalues: True
   
   genetic_algorithm:
     instance_count: 10
     generation_count: 50
     generation_size: 1000
     random_object_fraction: 0.2
     mutation_probability: 0.9
     crossover_probability: 0.25
     fitness_coefficient: 5.0
     verbose: True
   ```

   

4. leave R environment and run the following command in shell

   ```shell
   $ runSchism prepare_for_hypothesis_test -c W.yaml
   ```

   input:

   - `W.yaml`
   - `W.clusterEstimates.tsv`
   - ` W.mutation-to-cluster.tsv`

   output:

   - `0_W_results.cluster.cellularity`(needed)

5. then we use the Hypothesis test function in Schism

   ```shell
   $ runSchism hypothesis_test -c W.yaml
   ```

   input:

   - `0_W_results.cluster.cellularity`

   output:

   - `0_W_results.HT.cpov`( possibly needed)
   - `0_W_results.HT.pvalues`

6. and then the following command will give the evolution result of clusters

   ```shell
   $ runSchism run_ga --config W.yaml --mode serial
   $ runSchism summarize_ga_results -c W.yaml
   $ runSchism consensus_tree -c W.yaml
   ```

   input:

   - `0_W_results.HT.cpov`
   - `0_W_results.HT.pvalues`

   output:

   - `0_W_results.GA.consensusTree`(needed)

7. check  the directory of `0_W_results.GA.consensusTree` and `0_W_results.cluster.cellularity`

   ```
   dir.GA.consensusTree
   dir.cluster.cellularity
   ```

8. run the `schism2Fishplot` function in the same R environment and finish your fish-plotting

   ```R
   schism2Fishplot(dir.cluster.cellularity, dir.GA.consensusTree)
   ```