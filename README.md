<!--
 * @Autor: caoyh
 * @Date: 2022-11-17 09:27:27
 * @LastEditTime: 2022-11-18 10:10:54
-->
# MooSeeker

## Introduction
MooSeeker is a metabolic pathway design tool based on the multi-objective optimization algorithm that aims to trade off all the criteria optimally. 

- The metabolic pathway design problem is characterized as a multi-objective optimization problem with three objectives including pathway length, thermodynamic feasibility and theoretical yield. 
- In order to digitize the continuous metabolic pathway, MooSeeker develops the encoding strategy, BioCrossover and BioMutation operators to search for the candidate pathways. 
- Finally, MooSeeker outputs the Pareto optimal solutions of the candidate metabolic pathways with three criterion values.

The experiment results show that MooSeeker is capable of constructing the experimentally validated pathways and finding the higher-performance pathway than the single-objective-based methods. 

![The overveiw of MooSeeker](images/overall.jpg)

---

## Methodology
MooSeeker is proposed for design metabolic pathway with 3 stages:
### 1. Initilize Pairs Pool based on the KEGG Database 

```bash
python utils/init_pool.py
```
### 2. Improved NSGA-II algorithm for Metabolic Pathway 
The improved NSGA-II algorithm based on the **encoding strategy**, **BioCrossover** and **BioMutation** operators is shown in `model/My_NSGAII.py`.

### 3. Train the model 
According to different experiment, there are four types of config file for `multi-glycolysis.yaml`, `multi-vanillin.yaml`, `single-glycolysis.yaml`, `single-vanillin.yaml`.
We can train the model based on different config file for different result.

For example:
```bash
python main.py --config config/multi-vanillin.yaml
```


## Reference

[1] Yahui Cao, Tao Zhang, Xin Zhao, el. MooSeeker
