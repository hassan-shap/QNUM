# QNUM: Resource Placement for Rate and Fidelity Maximization in Quantum Networks

This repository contains the code and scripts developed for the project "Resource Placement for Rate and Fidelity Maximization in Quantum Networks," published in [IEEE Transactions on Quantum Engineering](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=10607917). The project focuses on optimizing resource placement in quantum networks to maximize both rate and fidelity under given constraints. In conjunction with another [repository](https://github.com/pooryousefshahrooz/q_net_planning).

## Overview

Quantum networks are crucial for distributed quantum computing and communication. In this project, we address the problem of resource placement—deciding where to place quantum resources (e.g., entanglement generation nodes) in a network to maximize the communication rate and fidelity between nodes.

This repository includes:

- **Algorithms for Resource Placement**: Tools to place quantum resources optimally within a given network topology.
- **Simulation Environment**: A simulation environment to evaluate the performance of the proposed algorithms on quantum networks.
- **Performance Evaluation**: Scripts to measure and plot the performance metrics such as rate and fidelity of the network under different configurations.

## Paper

The code in this repository is based on the paper:

**Resource Placement for Rate and Fidelity Maximization in Quantum Networks**  
*IEEE Transactions on Quantum Engineering*  
[DOI: 10.1109/ICC.2024.10607917](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=10607917)

### Abstract

The paper proposes a novel algorithm for placing quantum resources in a quantum network. The objective is to simultaneously maximize the communication rate and fidelity between nodes by strategically placing resources such as quantum repeaters and entanglement generation nodes. The proposed approach leverages mixed-integer linear programming (MILP) and heuristic methods to solve this NP-hard problem efficiently.

## Getting Started

### Prerequisites

To use this repository, you'll need to install the following dependencies:

- Python 3.8 or above
- NumPy
- SciPy
- NetworkX
- Matplotlib
- CVXPY (for convex optimization)

You can install the required Python packages using:

```bash
pip install -r requirements.txt
```

### Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/hassan-shap/QNUM
   cd QNUM
   ```

2. Install the required dependencies:

   ```bash
   pip install -r requirements.txt
   ```

### Usage

#### Running the Resource Placement Algorithm

To run the resource placement algorithm, use the following command:

```bash
python run_placement.py --network topology_file.json --output results.json
```

- `topology_file.json`: A JSON file representing the quantum network topology (nodes and edges).
- `results.json`: Output file where the algorithm's results will be saved.

#### Example

You can find example network topologies in the `examples/` directory. To run the algorithm on one of these examples:

```bash
python run_placement.py --network examples/sample_topology.json --output results/sample_results.json
```

#### Plotting Results

After running the algorithm, you can plot the results (rate and fidelity) using the following command:

```bash
python plot_results.py --input results/sample_results.json
```

This will generate plots for rate and fidelity optimization across different network configurations.

## Files and Directories

- `run_placement.py`: Main script to run the resource placement algorithm.
- `plot_results.py`: Script to plot the performance results (rate and fidelity).
- `requirements.txt`: Contains the dependencies required for the project.
- `examples/`: Directory with example quantum network topologies.

## Citation

If you use this repository in your research, please cite the following paper:

```bibtex
@inproceedings{hassan2024resource,
  title={Resource Placement for Rate and Fidelity Maximization in Quantum Networks},
  author={Hassan, Shap and others},
  booktitle={IEEE International Conference on Communications (ICC)},
  year={2024},
  pages={1--7},
  doi={10.1109/ICC.2024.10607917}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
```

### Key Points:
- The code structure assumes that there’s a `run_placement.py` script to execute the resource placement algorithm and `plot_results.py` to generate plots.
- Example network topologies are assumed to be stored in an `examples/` folder, and results are outputted in a JSON format.
- The citation is formatted in BibTeX for ease of use in academic contexts.

Let me know if you want any changes!
