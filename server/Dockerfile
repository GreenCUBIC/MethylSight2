FROM node:24-alpine3.21

EXPOSE 8989
WORKDIR /app
COPY . /app

RUN npm install

CMD ["npm", "run", "dev"]