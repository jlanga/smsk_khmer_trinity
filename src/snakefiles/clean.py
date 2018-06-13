rule clean:
    """
    Delete everything
    """
    shell:
        """
        if [ -d results ]; then
            rm -r results
        fi
        """
