# Work directory #
find /wrk/$USER -type d -not -path ".git/" -print0 | xargs -0 chmod u=rwx,g=rx,o=
find /wrk/$USER -type f -not -path ".git/" -print0 | xargs -0 chmod u=rw,g=r,o=

# Databases #
find ~/databases -type d -print0 | xargs -0 chmod u=rwx,g=rx,o=
find ~/databases -type f -print0 | xargs -0 chmod u=rw,g=r,o=