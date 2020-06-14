# import pymol
import numpy as np
from sys import argv
import argparse
import copy
import json


def calc_softmax_probs(x):
    """Calculate softmax probabilities from logits."""
    x_exp = np.exp(x - np.max(x, axis=-1, keepdims=True))
    return x_exp / np.sum(x_exp, axis=-1, keepdims=True)


def make_zero_diag(x):
    return x - np.diag(np.diag(x))


def pdb2map(pdb_filename):
	coords = []
	with open(pdb_filename) as pdbfile:
		for line in pdbfile:
			if (line[:4] == 'ATOM' or line[:6] == "HETATM") and line[12:16].strip() == "CA":
				x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
				coords.append(np.array([x, y, z]))
	seqlen = len(coords)
	contact_map = np.ones((seqlen, seqlen)) * -1
	for i, atom_a_xyz in enumerate(coords):
		for j, atom_b_xyz in enumerate(coords):
			distance = np.linalg.norm(atom_a_xyz - atom_b_xyz)
			contact_map[i, j] = distance
			contact_map[j, i] = distance
	return contact_map


def get_CA_sequence(ca_pdb_file):
	symbol_map = {
	    "ALA": "A",
	    "ARG": "R",
	    "ASN": "N",
	    "ASP": "D",
	    "CYS": "C",
	    "GLN": "Q",
	    "GLU": "E",
	    "GLY": "G",
	    "HIS": "H",
	    "ILE": "I",
	    "LEU": "L",
	    "LYS": "K",
	    "MET": "M",
	    "PHE": "F",
	    "PRO": "P",
	    "SER": "S",
	    "THR": "T",
	    "TRP": "W",
	    "TYR": "Y",
	    "VAL": "V"
	}
	sequence = ""
	with open(ca_pdb_file, 'r') as file:
		for line in file:
			if (line[:4] == 'ATOM' or line[:6] == "HETATM") and line[12:16].strip() == "CA":
				if len(sequence) == 0:
					start_idx = int(line[23:26])
				sequence += symbol_map[line[17:20]]
	return sequence, start_idx


def npz2map(filename):
	with np.load(filename) as predictions:
		if "dcls-ca" in predictions:
			preds = predictions["dcls-ca"]
		elif "dist-ca" in predictions:
			preds = predictions["dist-ca"]

	nb_bins = 36
	dist_val_min = 2.0
	dist_val_max = 20.0
	bin_size = (dist_val_max - dist_val_min) / nb_bins
	dist_vals = dist_val_min + bin_size / 2 + bin_size * np.arange(nb_bins)

	seq_len = preds.shape[0]
	probs = calc_softmax_probs(preds)
	probs = (probs + np.transpose(probs, [1, 0, 2])) / 2.0  # make it symmetric
	probs_cntc = probs[:, :, 1:nb_bins + 1]  # contacting bins
	dist_mat = np.sum(dist_vals[None, None, :] * probs_cntc, axis=-1) / (np.sum(probs_cntc, axis=-1) + 1e-8)
	dist_mat = make_zero_diag(dist_mat)

	return dist_mat

def get_deviations(crystal_pdb_f, pred_npz_f, ca_pdb_f, threshold):
	sequence, start_idx = get_CA_sequence(crystal_pdb_f)

	# Get all distance matrices
	crystal_mat = pdb2map(crystal_pdb_f)
	dist_mat = npz2map(pred_npz_f)
	ca_mat = pdb2map(ca_pdb_f)

	full_length = dist_mat.shape[0]

	# Crystal matrix is usually offset by a few residues so we need to adjust for that
	if dist_mat.shape != crystal_mat.shape:
		padded_crystal_mat = np.zeros_like(dist_mat)
		start = start_idx
		end = start + crystal_mat.shape[0]
		padded_crystal_mat[start:end, start:end] = crystal_mat
	else:
		start = 0
		padded_crystal_mat = crystal_mat

	# Add a mask to zero out residues not in the crystal PDB
	mask = (padded_crystal_mat > 0).astype(np.float32)
	# Distance predictions seem to only predict until 20A
	dist_mask = (padded_crystal_mat < 20).astype(np.float32)
	
	# Calculate deviations, basically we take mean absolute deviations along rows (accounting for zeros)
	# crystal_dist_deviations gets deviations between crystal PDB and predicted distances
	# ca_dist_deviations gets deviations between CA trace PDB and predicted distances
	crystal_dist_deviations = np.abs(dist_mat - padded_crystal_mat) * mask * dist_mask
	crystal_dist_deviations = np.where(crystal_dist_deviations > threshold)
	ca_dist_deviations = np.abs(dist_mat - ca_mat) * mask * dist_mask
	ca_dist_deviations = np.where(ca_dist_deviations > threshold)
	
	return start, crystal_dist_deviations, ca_dist_deviations

def visualize(prot_id, threshold):
	crystal_pdb_file = "./files/{}.pdb".format(prot_id)
	pred_npz_file = "./files/{}.npz".format(prot_id)
	if pred_npz_file not in os.listdir("files"):
		pred_npz_file = "./files/{}_tf_std.npz".format(prot_id)
	if pred_npz_file not in os.listdir("files"):
		print("NPZ file should either by <prot_id>.npz or <prot_id>_tf_std.npz")
	ca_pdb_file = "./files/{}_ca.pdb".format(prot_id)

	sequence, start_idx = get_CA_sequence(crystal_pdb_file)
	full_sequence, _ = get_CA_sequence(ca_pdb_file)
	full_length = len(full_sequence)

	start, crystal_dist_deviations, ca_dist_deviations = get_deviations(crystal_pdb_file, pred_npz_file, ca_pdb_file, threshold)
	
	# Initialize names
	crystal_name = prot_id[-5:]
	ca_name = prot_id[-5:] + "_catrace"
	crystal_viol_name = prot_id[-5:] + "_crystal-dist_violations"
	ca_viol_on_crystal_name = prot_id[-5:] + "_catrace-dist_violations_on_crystal"
	ca_viol_name = prot_id[-5:] + "_catrace-dist_violations"

	# Crystal PDB (white)
	cmd.load(crystal_pdb_file, crystal_name)
	cmd.color("white", crystal_name)
	cmd.show("ribbon", crystal_name)
	
	# Crystal vs Distance Predictions (magenta)
	cmd.load(crystal_pdb_file, crystal_viol_name)
	cmd.color("magenta", crystal_viol_name)
	cmd.hide("everything", crystal_viol_name)
	drawn = []
	for i, x in enumerate(crystal_dist_deviations[0]):
		if (x, crystal_dist_deviations[1][i]) not in drawn:
			cmd.bond("resi {} and name CA and {}".format(x, crystal_viol_name), "resi {} and name CA and {}".format(crystal_dist_deviations[1][i], crystal_viol_name))
			cmd.show("sticks", "(resi {}, resi {}) and name CA and {}".format(x, crystal_dist_deviations[1][i], crystal_viol_name))
			drawn.append((x, crystal_dist_deviations[1][i]))
			drawn.append((crystal_dist_deviations[1][i], x))
	cmd.set("stick_color", "magenta", crystal_viol_name)
	cmd.set("stick_transparency", 0.5, crystal_viol_name)

	# CA Trace vs Distance Predictions, superimposed onto Crystal Structure (red)
	cmd.load(crystal_pdb_file, ca_viol_on_crystal_name)
	cmd.color("red", ca_viol_on_crystal_name)
	cmd.hide("everything", ca_viol_on_crystal_name)
	drawn = []
	for i, x in enumerate(ca_dist_deviations[0]):
		if (x, ca_dist_deviations[1][i]) not in drawn:
			cmd.bond("resi {} and name CA and {}".format(x, ca_viol_on_crystal_name), "resi {} and name CA and {}".format(ca_dist_deviations[1][i], ca_viol_on_crystal_name))
			cmd.show("sticks", "(resi {}, resi {}) and name CA and {}".format(x, ca_dist_deviations[1][i], ca_viol_on_crystal_name))
			drawn.append((x, ca_dist_deviations[1][i]))
			drawn.append((ca_dist_deviations[1][i], x))
	cmd.set("stick_color", "red", ca_viol_on_crystal_name)
	cmd.set("stick_transparency", 0.5, ca_viol_on_crystal_name)

	# CA Trace (yellow)
	cmd.load(ca_pdb_file, ca_name)
	cmd.align(ca_name, crystal_name)
	cmd.color("yellow", ca_name)
	cmd.show("ribbon", ca_name)

	# CA Trace vs Distance Predictions (cyan)
	cmd.load(ca_pdb_file, ca_viol_name)
	cmd.align(ca_viol_name, crystal_name)
	cmd.color("cyan", ca_viol_name)
	cmd.hide("everything", ca_viol_name)
	drawn = []
	for i, x in enumerate(ca_dist_deviations[0]):
		if (x, ca_dist_deviations[1][i]) not in drawn:
			cmd.bond("resi {} and name CA and {}".format(x, ca_viol_name), "resi {} and name CA and {}".format(ca_dist_deviations[1][i], ca_viol_name))
			cmd.show("sticks", "(resi {}, resi {}) and name CA and {}".format(x, ca_dist_deviations[1][i], ca_viol_name))
			drawn.append((x, ca_dist_deviations[1][i]))
			drawn.append((ca_dist_deviations[1][i], x))
	cmd.set("stick_color", "cyan", ca_viol_name)
	cmd.set("stick_transparency", 0.5, ca_viol_name)

	cmd.center(crystal_name)

	print("Done.")


# parse input arguments
parser = argparse.ArgumentParser(description='Visualize violations')
parser.add_argument('-p', '--prot_id', type=str, required=True,
                    help='Name of the PDB file, without the extension')
parser.add_argument('-t', '--threshold', default=10, type=float,
                    help='Violation threshold, in angstroms')
args = parser.parse_args()

use_angle_constraints = True
filter_violations = True

visualize(args.prot_id, args.threshold)
