USE test_iter1;

SELECT qview.name as query, rview.name as reference, GenomeSimilarity.similarity
	FROM GenomeSimilarity
    
	# double join the asm table
	JOIN Assemblies as qasm
		ON GenomeSimilarity.query = qasm.id
	JOIN Assemblies as rasm
		ON GenomeSimilarity.reference = rasm.id
        
	# double join the DanglView
	JOIN DanglView as qview
		ON qasm.main = qview.id
	JOIN DanglView as rview
		ON rasm.main = rview.id
