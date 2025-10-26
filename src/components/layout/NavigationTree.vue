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
import categoryData from '@/content/categories.json';

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
      // 根据当前路由，按 categories.json 原顺序构建导航树
      // currentPath 形如 notes/Programming/python/biopython/biopython
      const path = this.currentPath;
      if (!path) { this.navigationTree = []; return; }

      const segs = path.split('/').filter(Boolean);
      if (segs.length < 2) { this.navigationTree = []; return; }

      const type = segs[0];  // notes | projects | topics
      const group = segs[1]; // Omics | Programming | biocrawler | demo ...

      // 1) 在 categories.json 中定位到对应 item（保持 JSON 顺序）
      let targetItem = null;
      if (Array.isArray(categoryData)) {
        outer:
        for (const section of categoryData) {
          if (!Array.isArray(section.items)) continue;
          for (const item of section.items) {
            if (item?.name !== group) continue;

            // 确认该 item 下存在属于当前 type 的文章（兼容老结构与新结构）
            let hasTypeMatch = false;

            if (Array.isArray(item.articles)) {
              hasTypeMatch = item.articles.some(a => typeof a?.articleUrl === 'string' && a.articleUrl.includes(`/article/${type}/`));
            }
            if (!hasTypeMatch && Array.isArray(item.categories)) {
              hasTypeMatch = item.categories.some(cat =>
                Array.isArray(cat.articles) &&
                cat.articles.some(a => typeof a?.articleUrl === 'string' && a.articleUrl.includes(`/article/${type}/${group}/`))
              );
            }

            if (hasTypeMatch) {
              targetItem = item;
              break outer;
            }
          }
        }
      }
      if (!targetItem) { this.navigationTree = []; return; }

      // 2) 生成根层文件（老结构）与子分类（新结构），严格按 JSON 顺序
      const rootFiles = [];
      const children = [];

      const toFile = (title, articleUrl) => {
        const parts = String(articleUrl).replace(/^\/+/, '').split('/');
        const i0 = parts[0] === 'article' ? 1 : 0;
        const t = parts[i0];
        const g = parts[i0 + 1];
        const rest = parts.slice(i0 + 2); // [subKey, ..., fileName]
        const pathNoExt = `${t}/${g}/${rest.join('/')}`;
        return { title, path: `${pathNoExt}.md` };
      };

      // 老结构：item.articles 直接挂根层
      if (Array.isArray(targetItem.articles)) {
        targetItem.articles.forEach(a => { if (a?.articleUrl) rootFiles.push(toFile(a.title, a.articleUrl)); });
      }

      // 新结构：item.categories[].articles -> children
      if (Array.isArray(targetItem.categories)) {
        targetItem.categories.forEach(cat => {
          const files = [];
          if (Array.isArray(cat.articles)) {
            cat.articles.forEach(a => { if (a?.articleUrl) files.push(toFile(a.title, a.articleUrl)); });
          }
          // 只在该子分类下有文件时渲染
          if (files.length) {
            children.push({
              name: cat.title || cat.key,
              type: 'directory',
              files
            });
          }
        });
      }

      // 3) 设定导航树根节点。根名使用 item.title（如“组学”、“编程”等）
      this.navigationTree = [{
        name: targetItem.title || targetItem.name || group,
        files: rootFiles,
        children
      }];
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
  color: var(--app-text-muted);      /* 更浅灰色 */
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
  color: var(--app-text-emphasis);
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
  color: var(--app-text-muted);
  border-radius: 0.25rem;
  transition: background-color 0.2s ease, color 0.2s ease;
}
.article-link:hover,
.article-link:focus {
  background-color: var(--app-primary-bg-subtle);
  color: var(--app-primary);
}
.article-link.active {
  background-color: transparent;  /* 去背景，强调纯文字样式 */
  color: var(--app-primary);       /* 蓝色 */
  font-weight: 700;               /* 粗体 */
}
</style>
