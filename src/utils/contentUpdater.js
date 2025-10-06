const fs = require('fs');
const path = require('path');
const matter = require('gray-matter');

const contentDir = path.join(__dirname, '../content');
const notesDir = path.join(contentDir, 'notes');

function getFileMetadata(filePath) {
    const fileContent = fs.readFileSync(filePath, 'utf-8');
    const { data } = matter(fileContent);
    const fileName = path.basename(filePath);
    const title = data.title || fileName.replace(/\.md$/, '').replace(/^\d{2}-\d{2}-\d{2}--/, '');
    return { title, path: path.relative(contentDir, filePath).replace(/\\/g, '/') };
}

function buildDirectoryTree(dir) {
    const items = fs.readdirSync(dir, { withFileTypes: true });
    const children = [];
    const files = [];

    for (const item of items) {
        const fullPath = path.join(dir, item.name);
        if (item.isDirectory()) {
            children.push(buildDirectoryTree(fullPath));
        } else if (item.isFile() && item.name.endsWith('.md')) {
            files.push(getFileMetadata(fullPath));
        }
    }

    const result = { name: path.basename(dir), type: 'directory' };
    if (children.length) result.children = children;
    if (files.length) result.files = files;
    return result;
}

function updateNotes() {
    const notesTree = buildDirectoryTree(notesDir).children;
    const outputPath = path.join(notesDir, 'notes.json');
    fs.writeFileSync(outputPath, JSON.stringify({ notes: notesTree }, null, 2));
    console.log('notes.json has been updated successfully.');
}

function updateProjects() {
    // Logic for updating projects can be added here
}

function updateTopics() {
    // Logic for updating topics can be added here
}

function updatePosts() {
    // Logic for updating posts can be added here
}

// Run the update
updateNotes();