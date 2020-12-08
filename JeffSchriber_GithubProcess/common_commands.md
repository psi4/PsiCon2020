## Common git commands:
  0. `git clone <https://github.com/some_repo/some_repo.git>` Obtain exact copy of repo from github
  1. `git add <file>` Adds local file to git tracking
  2. `git log` Shows local history of commits
  3. `git status` Shows any changes to local files, in addition to any untracked files
  4. `git commit -m "text message"` Commits any files that have been added to stage with message
  5. `git commit -a` Adds then commits all changes in tracked files
  6. `git pull <remote> <branch>` Merges github history to local, taken from remote (usually origin) and branch (often master).
  7. `git push <remote> <branch>` Moves local commits to github repository at remote/branch
  8. `git checkout -b <branch_name>` Makes a new branch and switch to it
  9. `git checkout <branch_name>` Move to branch_name
  10. `git branch` Lists branches, shows which branch you are on
  11. `git branch -d` Delete branch locally
  12. `git push <remote> :<branch_name>` Deletes branch on remote.
  13. `git checkout <file>` Reverts any changes in `file` to state determined by current git hash
  14. `git init` Create a new git repository
  15. `git remote add origin <server>` Connects repository to a remote server
  16. `git diff` View all local changes relative to remote
  17. `git rm` Removes a tracked file
  18. `git remote -v` Shows identity of origin and other remotes
  19. `git remote add upstream <parent repo>` Adds upstream remote for a parent repo


## Making a new repo
1. First make repo on github
2. `git clone` your new repo locally.
3. Make some new file, `file.dat`
4. See your new untracked file: `git status`
5. Stage file for commit: `git add file.dat`
6. Commit change, enter message, close file: `git commit -a`
7. See your new history!: `git log`
8. Push your changes to your repo on github.com: `git push origin master`

## Making a fork
### Good practice for making changes to shared code
1. On a parent repo page in github.com, click Fork on the top right.
2. Clone to local machine: `git clone <your fork>`
3. Setup upstream: `git remote add upstream <parent repository>`
4. Make a branch of your fork: `git checkout -b new_branch` (branch name should be descriptive of changes)
5. Make some local change: `vi new_file.dat`
6. Add and commit change: `git add new_file.dat`, `git commit -m "Add new file"`
7. Rebase to resolve any merge conflicts: `git pull --rebase upstream master`
8. If your local branch differs from origin, you will need to make your changes and then push with force:
`git push --force origin new_branch`
9. Push changes: `git push origin new_branch` (this will create a new branch on github.com)
10. On your fork on github.com, create a pull request pushing the "Compare & pull request" button.
11. Open the pull request by by writing a description and clicking "Create pull request"
12. Delete branch once PR is merged
