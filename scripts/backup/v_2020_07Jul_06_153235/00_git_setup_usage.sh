# Fork repository from the github
# https://github.com/nf-core/atacseq

# Set config values
git config user.name  "Gaurav Jain"
git config user.email "gauravj49@gmail.com"

# Check the config list
git config --list

# Get the status of the project and repository
git status

# Ignore files that should not go into the repository
# emacs -nw .gitignore
cat > .gitignore
annotation
docs
input
output
.gitignore
00_git_setup_usage.sh

# Once you have done that, git now knows about your remote repository. You can then tell it to push (which is "upload") your commited files:
git push
