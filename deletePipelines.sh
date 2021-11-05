#!/bin/bash
set -e

GITLAB_URL="https://gitlab.dlr.de"
TOKEN=""
PROJECT="13683"
# How many to delete from the oldest.
PER_PAGE=100
UPDATED_BEFORE=2021-07-01T00:00:00Z

echo "The following pipelines will be deleted: "
for PIPELINE in $(curl -s --header "PRIVATE-TOKEN: $TOKEN" "$GITLAB_URL/api/v4/projects/$PROJECT/pipelines?per_page=$PER_PAGE&sort=asc&updated_before=$UPDATED_BEFORE" | jq '.[].id') ; do
    echo "$PIPELINE"
done

read -p "Are you sure? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    for PIPELINE in $(curl -s --header "PRIVATE-TOKEN: $TOKEN" "$GITLAB_URL/api/v4/projects/$PROJECT/pipelines?per_page=$PER_PAGE&sort=asc&updated_before=$UPDATED_BEFORE" | jq '.[].id') ; do
        echo "Deleting pipeline $PIPELINE"
        curl --header "PRIVATE-TOKEN: $TOKEN" --request "DELETE" "https://gitlab.dlr.de/api/v4/projects/$PROJECT/pipelines/$PIPELINE"
    done
fi
