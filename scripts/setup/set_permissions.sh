find ~/work/ -type d -print0 | xargs -0 chmod u=rwx,g=rx,o=
find ~/work/ -type f -print0 | xargs -0 chmod u=rw,g=r,o=
