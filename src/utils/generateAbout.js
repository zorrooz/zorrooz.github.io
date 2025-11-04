import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import yaml from 'js-yaml';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const contentSrcDir = path.join(__dirname, '../content-src');
const contentOutputDir = path.join(__dirname, '../content');

function ensureDirectoryExistence(filePath) {
  const dirname = path.dirname(filePath);
  if (!fs.existsSync(dirname)) {
    fs.mkdirSync(dirname, { recursive: true });
  }
}

function normalize(raw = {}) {
  const intro = typeof raw.introduction === 'string' ? raw.introduction : '';

  const section = Array.isArray(raw.section)
    ? raw.section.map(s => {
        const title = typeof s?.title === 'string' ? s.title : '';
        const items = Array.isArray(s?.items)
          ? s.items.map(it => ({
              item: typeof it?.item === 'string' ? it.item : '',
              desc: typeof it?.desc === 'string' ? it.desc : ''
            }))
          : [];
        return { title, items };
      })
    : [];

  const contacts = Array.isArray(raw.contacts) ? raw.contacts : [];

  return {
    introduction: intro,
    section,
    contacts
  };
}

function generateAboutJson(locale = 'zh-CN') {
  const yamlPath = path.join(contentSrcDir, locale === 'zh-CN' ? 'about.yaml' : 'about-en.yaml');
  const outputPath = path.join(contentOutputDir, locale === 'zh-CN' ? 'about.json' : 'about-en.json');

  let raw = {};
  if (fs.existsSync(yamlPath)) {
    try {
      const yamlContent = fs.readFileSync(yamlPath, 'utf-8');
      raw = yaml.load(yamlContent) || {};
    } catch (e) {
      console.warn(`Warn: failed to parse YAML at ${yamlPath}, using empty object.`, e?.message || e);
      raw = {};
    }
  } else {
    console.warn(`Warn: source YAML not found at ${yamlPath}, using empty object.`);
  }

  const result = normalize(raw);

  try {
    ensureDirectoryExistence(outputPath);
    fs.writeFileSync(outputPath, JSON.stringify(result, null, 2), 'utf-8');
    console.log(`Successfully generated: ${outputPath}`);
  } catch (error) {
    console.error(`Failed to generate ${outputPath}:`, error);
  }
}

function main() {
  console.log('Starting about.json generation script...');
  generateAboutJson('zh-CN');
  generateAboutJson('en-US');
  console.log('about.json generation complete.');
}

main();

export { generateAboutJson };