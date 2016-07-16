# PostgreSQL

## How to connect

`/etc/default/portal-frontend`

```bash
psql postgres://binatech:binarocks2013@localhost/portal
```
## List all databases

```bash
=> \l
```

## Quit psql
```bash
=> \q
```

## How to drop database `portal`

```bash
sudo systemctl stop portal-frontend&& sudo systemctl stop portal-backend && sudo systemctl stop executor

sudo -u postgres bash

dropdb portal
```

Recreate the database

```bash
sudo salt-call apply postgres

sudo systemctl restart portal-frontend && sudo systemctl restart portal-backend && sudo systemctl restart executor
```

