goosehints:
  [ -f .goosehints ] || ln -s CLAUDE.md .goosehints

copilot-instructions:
  [ -f .github/copilot-instructions.md ] || cd .github && ln -s ../CLAUDE.md copilot-instructions.md

setup-ai: setup-ai-instructions setup-gh

setup-ai-instructions: goosehints copilot-instructions

setup-gh: gh-add-topics gh-add-secrets

gh-add-topics:
  gh repo edit --add-topic "monarchinitiative,ai4curation"

gh-add-secrets:
  gh secret set PAT_FOR_PR --body "$PAT_FOR_PR"
  gh secret set ANTHROPIC_API_KEY --body "$ANTHROPIC_API_KEY"
  gh secret set OPENAI_API_KEY --body "$OPENAI_API_KEY"
  gh secret set CBORG_API_KEY --body "$CBORG_API_KEY"
  gh secret set CLAUDE_CODE_OAUTH_TOKEN --body "$GH_CLAUDE_CODE_OAUTH_TOKEN"

