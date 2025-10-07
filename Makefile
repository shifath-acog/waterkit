JP_IMAGE := waterkit-cli
JP_CONTAINER := waterkit-cli
JP_PORT := $(shell shuf -i 8000-17999 -n 1)

SERVICE_NAME := waterkit

build-image:
	@cd devops && \
		JP_IMAGE=$(JP_IMAGE) JP_CONTAINER=$(JP_CONTAINER) JP_PORT=$(JP_PORT) \
		docker compose build --no-cache

start-container:
	@docker ps --format '{{.Names}}' | grep -q "^$(JP_CONTAINER)$$" && echo "Already running container: \e[1;32m$(JP_CONTAINER)\e[0m" || true

	@if ! docker ps --format '{{.Names}}' | grep -q -e "^$(JP_CONTAINER)$$"; then \
		cd devops && \
			JP_IMAGE=$(JP_IMAGE) JP_CONTAINER=$(JP_CONTAINER) JP_PORT=$(JP_PORT) \
			docker compose -p $(SERVICE_NAME) up && \
		echo "Successfully completed simulation: \e[1;32m$(JP_CONTAINER)\e[0m"; \
	fi

stop-container:
	@if docker ps -a --format '{{.Names}}' | grep -q -e "^$(JP_CONTAINER)$$"; then \
		echo "Stopped and removed container: \033[1;31m$(JP_CONTAINER)\033[0m"; \
		docker rm -f $(JP_CONTAINER) > /dev/null 2>&1; \
	else echo "\033[1;31mThere are no running containers to stop.\033[0m"; \
	fi
