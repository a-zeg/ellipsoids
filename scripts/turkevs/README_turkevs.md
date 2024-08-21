This README contains instructions for running experiments analogous to https://arxiv.org/abs/2206.10551, including classification based on ellipsoid and Rips complexes.

These scripts make heavy use of the code available here [https://github.com/renata-turkes/turkevs2022on](https://github.com/renata-turkes/turkevs2022on).

To run the experiments, the following steps are necessary:
1. Generate the datasets by running the `generate_turkevs.py` script;
2. Calculate the ellipsoids (and Rips) barcodes by running the `calculate_turkevs.py` script;
3. Run the experiments by running the `run_turkevs_holes_experiments.py` script.