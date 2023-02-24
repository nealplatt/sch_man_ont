cd /master/nplatt/sm_nanopore

while true; do
        if test -z "$(git status -s)"; then
                echo "$(date): No changes found. Skipping push."
        else
                echo "$(date): Changes found. Pushing changes..."
                git add \*.py \*.ipynb \*.sh && git commit -m "auto-update $(date)" && git push
        fi

        sleep 24h
done
