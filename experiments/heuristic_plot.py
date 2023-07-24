import os
# import git  # python -m pip install gitpython
import matplotlib.pyplot as plt


def plot(filename, scores, nsh_scores):
    ax = plt.axes()
    ax.set_title(f'Comparison: {filename}')
    ax.set_ylabel('Best score')

    ax.plot(range(len(scores)), scores, label='Heuristic')
    ax.plot(range(len(nsh_scores)), nsh_scores, label='Non-heuristic')
    ax.legend()

    plt.show()


def main():
    # root = git.Repo('.', search_parent_directories=True).git.rev_parse("--show-toplevel")
    # dirs = [os.path.join(root, 'experiments', 'structure_recovery', 'results'),
            # os.path.join(root, 'experiments', 'real_world_graphs', 'results')]
    dirs = [os.path.join('.', 'structure_recovery', 'results'),
            os.path.join('.', 'real_world_graphs', 'results')]

    for rdir in dirs:
        files = [file for file in os.listdir(rdir)
                 if file.split('/')[-1][0:4] == 'nsh_'
                 and file.split('.')[0].split('_')[-1] not in ('graph', 'noise', 'nodes', 'runtimes')]

        for file in files:
            scores = []
            nsh_scores = []

            with open(os.path.join(rdir, file[4:]), 'r') as logfile, \
                 open(os.path.join(rdir, file), 'r') as nsh_logfile:

                for line in logfile:
                    if 'best score' in line.lower():
                        score = float(line.strip().split(' ')[-1])
                        scores.append(score)

                for line in nsh_logfile:
                    if 'best score' in line.lower():
                        nsh_score = float(line.strip().split(' ')[-1])
                        nsh_scores.append(nsh_score)

            if scores and nsh_scores:
                plot(file[4:], scores, nsh_scores)
            else:
                print(file)


if __name__ == '__main__':
    main()
