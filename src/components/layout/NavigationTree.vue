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
import notesFlat from '@/content/notes.json';

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
      // 从扁平 notes.json 动态构建 notesTree
      const buildTreeFromFlat = (items) => {
        const rootChildren = new Map(); // categoryName -> node
        const nameMap = {
          'Programming': '编程语言',
          'Bioinformatics': '生物信息学',
          'Omics': '组学技术',
          'DataScience': '数据科学',
          'python': 'Python',
          'r': 'R语言',
          'shell': 'Shell',
          'javascript': 'JavaScript',
          'alignment': '序列比对',
          'structure': '结构分析',
          'genomics': '基因组学',
          'proteomics': '蛋白质组学',
          'transcriptomics': '转录组学',
          'statistics': '统计分析',
          'machinelearning': '机器学习',
          'visualization': '数据可视化'
        };
        const formatName = (n) => nameMap[n] || n;

        const ensureNode = (parent, name) => {
          if (!parent.children) parent.children = [];
          let node = parent.children.find(n => n._rawName === name || n.name === formatName(name));
          if (!node) {
            node = { name: formatName(name), _rawName: name, type: 'directory', children: [], files: [] };
            parent.children.push(node);
          }
          return node;
        };

        const virtualRoot = { children: [] };
        if (Array.isArray(items)) {
          items.forEach(it => {
            const rel = it.relativePath || '';
            const segs = rel.split('/').filter(Boolean);
            if (segs.length === 0) return;

            // 如果末两段相同（如 bwa/bwa），跳过倒数第二层冗余目录
            const end = (segs.length >= 2 && segs[segs.length - 1] === segs[segs.length - 2])
              ? segs.length - 2
              : segs.length - 1;

            let current = ensureNode(virtualRoot, segs[0]);
            for (let i = 1; i <= end - 1; i++) {
              current = ensureNode(current, segs[i]);
            }

            // 文件节点
            const fileTitle = it.title || segs[segs.length - 1];
            const filePath = `notes/${rel}.md`;
            if (!current.files) current.files = [];
            current.files.push({ title: fileTitle, path: filePath });
          });
        }

        // 清理辅助字段并返回
        const clean = (node) => {
          if (!node) return node;
          const copy = { ...node };
          delete copy._rawName;
          if (copy.children && copy.children.length) {
            copy.children = copy.children.map(clean);
          } else {
            delete copy.children;
          }
          if (copy.files && copy.files.length === 0) delete copy.files;
          return copy;
        };
        return (virtualRoot.children || []).map(clean);
      };

      const notesTree = buildTreeFromFlat(notesFlat);
      const combinedTree = [...notesTree];

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
        const notesRootNodes = fullTree; // 仅 notes
        if (this.currentPath === 'notes') {
          this.navigationTree = notesRootNodes;
          return;
        }
        const activeRoot = findContainingRoot(notesRootNodes, this.currentPath);
        this.navigationTree = activeRoot ? [activeRoot] : [];
      } else {
        this.navigationTree = []; // 仅支持 notes
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
  padding-left: 0.75rem;

}

/* 根列表不缩进，让二级分类对齐一级分类 */
.article-list { list-style: none; padding-left: 0; }
.article-list-root { padding-left: 0; }

/* 二级分类：黑色正常体，并与一级对齐 */
.directory-node { padding: 0.25rem 0; }
.directory-name.level-2 {
  display: block;
  font-weight: 700;               /* 加粗 */
  color: var(--bs-body-color);    /* 黑色（正文颜色） */
  margin-left: 0;                 /* 与一级对齐 */
  padding-left: 0.75rem;
}

/* 三级文章：有缩进，黑色或灰色；当前文章蓝色粗体 */
.sub-list { padding-left: 0.75rem; margin-top: 0.5rem; }
.files-level { padding-left: 0.75rem; } /* 三级文章缩进 */

.article-item { margin-bottom: 0.3rem; }
.article-link {
  display: block;
  padding: 0.5rem 0.75rem;
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
