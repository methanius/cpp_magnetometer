import sys
import numpy as np
import pickle

PQS_av_widths = []
forward_av_widths = []
PQS_mean_errors = []
forward_mean_errors = []


for i in range(1, len(sys.argv)):
        with open(sys.argv[i], 'rb') as f:
                [weighted_t_averaged_width, mean_errorbar, weighted_t_averaged_width_forward, mean_errorbar_forward] = pickle.load(f)
                PQS_av_widths.append(weighted_t_averaged_width)
                forward_av_widths.append(weighted_t_averaged_width_forward)
                PQS_mean_errors.append(mean_errorbar)
                forward_mean_errors.append(mean_errorbar_forward)

PQS_av_widths = np.asarray(PQS_av_widths)
PQS_mean_errors = np.asarray(PQS_mean_errors)
forward_av_widths = np.asarray(forward_av_widths)
forward_mean_errors = np.asarray(forward_mean_errors)


final_PQS_width = PQS_av_widths.sum() / len(PQS_av_widths)
final_forward_width = forward_av_widths.sum() / len(PQS_av_widths)
final_PQS_error = PQS_mean_errors.sum() / len(PQS_mean_errors)
final_forward_error = forward_mean_errors.sum() / len(forward_mean_errors)

with open('averaged_widths.txt', 'w') as f:
    f.write('Averaged PQS width = {}\n'.format(final_PQS_width))
    f.write('Averaged PQS error = {}\n'.format(final_PQS_error))
    f.write('Averaged forward width = {}\n'.format(final_forward_width))
    f.write('Averaged forward error = {}\n'.format(final_forward_error))
