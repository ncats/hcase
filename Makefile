images:
	@echo "Building images..."
	@docker build -t hcase -f docker/Dockerfile.hcase .

jupyter:
	docker run -it --rm \
		-v $(PWD):/app \
		-p 8874:8874 hcase \
	    jupyter-lab --allow-root --port=8874 --ip 0.0.0.0 --no-browser --notebook-dir=/app



ruff:
	docker run -it --rm -v $(PWD):/app --env-file .env selective-inference-dev ruff .

isort:
	docker run -it --rm -v $(PWD):/app --env-file .env selective-inference-dev isort .

mypy:
	docker run -it --rm -v $(PWD):/app --env-file .env --workdir /app selective-inference-dev mypy src

pytest:
	docker run -it --rm -v $(PWD):/app --env-file .env selective-inference-dev ptw --runner python run_pytest.py .
