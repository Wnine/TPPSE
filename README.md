# TPPSE

Implementation of the **Three-Patterns-Protected Searchable Encryption Supporting Disjunctive Keyword Search(TPPSE)** scheme as described in our paper.

## Introduction

This repository provides an implementation of TPPSE (Three-Patterns-Protected Searchable Encryption).  
The scheme supports disjunctive keyword search and defends against multiple types of search pattern leakage, ensuring robust privacy protection.

## Leakage protection

The TPPSE scheme protects against:
- **SP (Search Pattern)** — Hides whether the same query is issued multiple times.
- **QLP (Query Length Pattern)** — Hides the number of keywords involved in each query (**n**).
- **AP (Access Pattern)** — Hides which documents match the query.
- **RLP (Response Length Pattern)** — Hides the number of documents returned in response to a query.

By concealing these patterns, the scheme prevents adversaries from inferring sensitive information about the queries and the dataset.

## Project structure

```
TPPSE/
├── Ashe_galois.py      # Finite field operations and polynomial encoding
├── search.py           # Protocol 2 search
├── setup.py            # Algorithm 1 setup
├── requirements.txt    # Python dependencies
├── LICENSE             # MIT License
├── .gitignore          # Git ignore configuration
```

## Installation

### Clone this repository
```bash
git clone https://github.com/Wnine/TPPSE.git
cd TPPSE
```

### Install dependencies
```bash
pip install -r requirements.txt
```

## Example usage

Run the main experimental script:
```bash
python search.py
```
*(Adjust as needed to provide appropriate input data or configurations. If database files are required, please contact the author to obtain them)*

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
