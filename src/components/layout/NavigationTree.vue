<template>
  <div class="navigation-tree">
    <div v-for="category in navigationTree" :key="category.name" class="category-group">
      <h3 class="category-name">{{ category.name }}</h3>
      <ul class="article-list article-list-root">
        <!-- 二级分类（与一级对齐，不缩进） -->
        <li v-if="category.children && category.children.length" v-for="dir in category.children" :key="dir.name" class="article-item">
          <div class="directory-node">
            <span class="directory-name level-2">{{ dir.name }}</span>
            <!-- 三级文章（有缩进） -->
            <ul v-if="dir.files && dir.files.length" class="article-list sub-list files-level">
              <li v-for="file in dir.files" :key="file.path" class="article-item">
                <router-link :to="toArticle(file.path)" class="article-link level-3" :class="{ active: isActive(file.path) }">
                  {{ file.title }}
                </router-link>
              </li>
            </ul>
            <!-- 更深层子目录 -->
            <ul v-if="dir.children && dir.children.length" class="article-list sub-list">
              <li v-for="sub in dir.children" :key="sub.name" class="article-item">
                <div class="directory-node">
                  <span class="directory-name level-2">{{ sub.name }}</span>
                  <ul v-if="sub.files && sub.files.length" class="article-list sub-list files-level">
                    <li v-for="file in sub.files" :key="file.path" class="article-item">
                      <router-link :to="toArticle(file.path)" class="article-link level-3" :class="{ active: isActive(file.path) }">
                        {{ file.title }}
                      </router-link>
                    </li>
                  </ul>
                </div>
              </li>
            </ul>
          </div>
        </li>
        <!-- 若一级分类下直接存在文件，也作为三级文章处理 -->
        <li v-if="category.files && category.files.length" v-for="file in category.files" :key="file.path" class="article-item">
          <router-link :to="toArticle(file.path)" class="article-link level-3" :class="{ active: isActive(file.path) }">
            {{ file.title }}
          </router-link>
        </li>
      </ul>
    </div>
  </div>
</template>

<script>
import notesData from '@/content/notes/notes.json';
import projectsData from '@/content/projects/projects.json';
import topicsData from '@/content/topics/topics.json';

export default {
  name: 'NavigationTree',
  data() {
    return {
      navigationTree: [],
      currentPath: ''
    };
  },
  watch: {
    '$route': {
      handler(to) {
        this.currentPath = to.params.path ? to.params.path.join('/') : '';
        this.buildTree();
      },
      immediate: true,
      deep: true
    }
  },
  mounted() {
    // buildTree is called by the route watcher on initialization
  },
  methods: {
    toArticle(path) {
      return { name: 'Article', params: { path: path.replace(/\.md$/, '').split('/') } };
    },
    isActive(path) {
      return this.currentPath === path.replace(/\.md$/, '');
    },
    buildTree() {
      // Build the full navigation data structure
      const notesTree = JSON.parse(JSON.stringify(notesData.notes));
      const projectsTree = { name: '项目', type: 'directory', children: JSON.parse(JSON.stringify(projectsData)) };
      const topicsTree = { name: '专题', type: 'directory', children: JSON.parse(JSON.stringify(topicsData)) };
      const combinedTree = [...notesTree, projectsTree, topicsTree];

      const processNode = (node) => {
        // 排序目录下的子目录与文件
        if (node.children) {
          node.children.forEach(processNode);
          node.children.sort((a, b) => (a.name || '').localeCompare(b.name || ''));
        }
        if (node.files) {
          node.files.sort((a, b) => (a.title || a.path).localeCompare(b.title || b.path));
        }
        return node;
      };
      const fullTree = combinedTree.map(processNode);

      // Filter the tree based on the current route
      if (!this.currentPath) {
        this.navigationTree = []; // Show nothing if not in a specific content path
        return;
      }

      const pathSegments = this.currentPath.split('/');
      const topLevel = pathSegments[0];

      const findContainingRoot = (nodes, targetPath) => {
        for (const node of nodes) {
          const searchDescendants = (currentNode) => {
            // 在 files 中查找
            if (currentNode.files && currentNode.files.length) {
              if (currentNode.files.some(f => f.path.replace(/\.md$/, '') === targetPath)) {
                return true;
              }
            }
            // 递归 children
            if (currentNode.children && currentNode.children.length) {
              return currentNode.children.some(child => searchDescendants(child));
            }
            return false;
          };
          if (searchDescendants(node)) {
            return node;
          }
        }
        return null;
      };

      if (topLevel === 'notes') {
        const notesRootNodes = fullTree.filter(node => node.name !== '项目' && node.name !== '专题');
        if (this.currentPath === 'notes') {
          this.navigationTree = notesRootNodes;
          return;
        }
        const activeRoot = findContainingRoot(notesRootNodes, this.currentPath);
        this.navigationTree = activeRoot ? [activeRoot] : [];
      } else if (topLevel === 'projects') {
        const projectsRoot = fullTree.find(node => node.name === '项目');
        if (projectsRoot) {
          const slug = pathSegments[1];
          const activeProject = projectsRoot.children.find(p => p.name === slug);
          if (activeProject) {
            this.navigationTree = [{ ...projectsRoot, children: [activeProject] }];
          } else {
            this.navigationTree = [projectsRoot];
          }
        } else {
          this.navigationTree = [];
        }
      } else if (topLevel === 'topics') {
        const topicsRoot = fullTree.find(node => node.name === '专题');
        if (topicsRoot) {
          const slug = pathSegments[1];
          const activeTopic = topicsRoot.children.find(t => t.name === slug);
          if (activeTopic) {
            this.navigationTree = [{ ...topicsRoot, children: [activeTopic] }];
          } else {
            this.navigationTree = [topicsRoot];
          }
        } else {
          this.navigationTree = [];
        }
      } else {
        this.navigationTree = []; // Default to showing nothing
      }
    }
  }
};
</script>

<style scoped>
.navigation-tree {
  padding: 0;
  font-size: 0.95rem;
}

/* 包裹每个一级分类组 */
.category-group { margin-bottom: 1rem; }

/* 一级分类：灰色粗体 */
.category-name {
  font-size: 1.1rem;
  font-weight: 700;               /* 粗体 */
  color: var(--bs-gray-600);      /* 更浅灰色 */
  margin-bottom: 0.5rem;
  padding-bottom: 0.25rem;

}

/* 根列表不缩进，让二级分类对齐一级分类 */
.article-list { list-style: none; padding-left: 0; }
.article-list-root { padding-left: 0; }

/* 二级分类：黑色正常体，并与一级对齐 */
.directory-node { padding: 0.25rem 0; }
.directory-name.level-2 {
  font-weight: 700;               /* 加粗 */
  color: var(--bs-body-color);    /* 黑色（正文颜色） */
  margin-left: 0;                 /* 与一级对齐 */
}

/* 三级文章：有缩进，黑色或灰色；当前文章蓝色粗体 */
.sub-list { padding-left: 0.75rem; margin-top: 0.5rem; }
.files-level { padding-left: 0.75rem; } /* 三级文章缩进 */

.article-item { margin-bottom: 0.3rem; }
.article-link {
  display: block;
  padding: 0.25rem 0.5rem;
  text-decoration: none;
  color: var(--bs-gray-700);      /* 默认灰/黑 */
  border-radius: 0.25rem;
  transition: background-color 0.2s ease, color 0.2s ease;
}
.article-link:hover {
  background-color: var(--bs-light);
  color: var(--bs-body-color);    /* 悬停黑色 */
}
.article-link.active {
  background-color: transparent;  /* 去背景，强调纯文字样式 */
  color: var(--bs-primary);       /* 蓝色 */
  font-weight: 700;               /* 粗体 */
}
</style>
